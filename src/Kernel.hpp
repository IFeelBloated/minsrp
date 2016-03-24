#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <vector>
#include <algorithm>
#include "VapourSynth.h"
#include "VSHelper.h"
using namespace std;

template<typename T>
static inline T _abs(T val) {
	return (val > 0) ? val : -val;
}

static inline void _clamp(long &val, int bits) {
	val = max(min(val, (1l << bits) - 1l), 0l);
}

template<typename T1, typename T2>
static void _pad(const T1 *&srcp, int &src_stride, T1 *&padp, int pad_stride, int h, int w) {
	T1 *padtmp = padp;
	for (int y = 1; y < h + 1; ++y) {
		padp += pad_stride;
		for (int x = 1; x < w + 1; ++x)
			padp[x] = srcp[x - 1];
		srcp += src_stride;
	}
	srcp = padtmp + pad_stride + 1;
	src_stride = pad_stride;
	padp += pad_stride;
	for (int x = 1; x < w + 1; ++x)
		padp[x] = padp[x - pad_stride];
	padp = padtmp;
	for (int x = 1; x < w + 1; ++x)
		padp[x] = padp[x + pad_stride];
	for (int y = 1; y < h + 1; ++y)
		padp[y*pad_stride] = padp[y*pad_stride + 1];
	padp = padp + pad_stride - 1;
	for (int y = 1; y < h + 1; ++y)
		padp[y*pad_stride] = padp[y*pad_stride - 1];
	padp = padtmp;
	*padp = static_cast<T1>((static_cast<T2>(padp[1]) + padp[pad_stride]) / 2);
	padp = padp + pad_stride - 1;
	*padp = static_cast<T1>((static_cast<T2>(*(padp - 1)) + padp[pad_stride]) / 2);
	padp += (h + 1)*pad_stride;
	*padp = static_cast<T1>((static_cast<T2>(*(padp - 1)) + *(padp - pad_stride)) / 2);
	padp = padtmp + (h + 1)*pad_stride;
	*padp = static_cast<T1>((static_cast<T2>(padp[1]) + *(padp - pad_stride)) / 2);
	padp = nullptr;
}

static void _gauss_3x3(float *&dstp, vector<float> &a, int x) {
	dstp[x] = static_cast<float>(
		(4. * a[0] +
			2. * (static_cast<double>(a[2]) + a[4] + a[5] + a[7]) +
			a[1] + a[3] + a[6] + a[8]) / 16.);
}

static void _gauss_3x3(uint16_t *&dstp, vector<uint16_t> &a, int x) {
	long val = ((a[0] << 2l) +
		((static_cast<uint32_t>(a[2]) + a[4] + a[5] + a[7]) << 1l) +
		a[1] + a[3] + a[6] + a[8] + 8l) >> 4l;
	_clamp(val, 16);
	dstp[x] = static_cast<uint16_t>(val);
}

static void _gauss_3x3(uint8_t *&dstp, vector<uint8_t> &a, int x) {
	long val = ((a[0] << 2l) +
		((static_cast<uint32_t>(a[2]) + a[4] + a[5] + a[7]) << 1l) +
		a[1] + a[3] + a[6] + a[8] + 8l) >> 4l;
	_clamp(val, 8);
	dstp[x] = static_cast<uint8_t>(val);
}

template<typename T>
static void _median_3x3(T *&dstp, vector<T> &a, int x) {
	nth_element(a.begin(), a.begin() + a.size() / 2, a.end());
	dstp[x] = a[a.size() / 2];
}

template<typename T>
static void _process_plane_3x3(const T *&srcp, int src_stride, T *&dstp, int dst_stride, int h, int w, bool gauss) {
	const T *srctmp = srcp;
	T *dsttmp = dstp;
	void(*_func)(T *&, vector<T> &, int) = nullptr;
	if (gauss)
		_func = _gauss_3x3;
	else
		_func = _median_3x3;
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			vector<T> a{ *(srcp + x),
				*(srcp + x - src_stride - 1), *(srcp + x - src_stride), *(srcp + x - src_stride + 1),
				*(srcp + x - 1), *(srcp + x + 1),
				*(srcp + x + src_stride - 1), *(srcp + x + src_stride), *(srcp + x + src_stride + 1) };
			_func(dstp, a, x);
		}
		dstp += dst_stride;
		srcp += src_stride;
	}
	srcp = srctmp;
	dstp = dsttmp;
}

static inline void _reverse_dif(const float *&srcp, float *&dstp, int x, double _amp, bool _linear) {
	if (_linear)
		dstp[x] = static_cast<float>((1. + _amp) * srcp[x] - dstp[x] * _amp);
	else
		dstp[x] = static_cast<float>(srcp[x] +
			(_amp * 4. * pow(_abs(srcp[x] * 256. - dstp[x] * 256.) / 4., 0.25)) *
			((srcp[x] * 256. - dstp[x] * 256.) / (_abs(srcp[x] * 256. - dstp[x] * 256.) + 1.001)) / 256.);
}

static inline void _reverse_dif(const uint16_t *&srcp, uint16_t *&dstp, int x, double _amp, bool _linear) {
	long val = 0l;
	if (_linear)
		val = static_cast<long>((1. + _amp) * srcp[x] - dstp[x] * _amp);
	else
		val = static_cast<long>(srcp[x] +
			(_amp * 4. * pow(_abs(srcp[x] / 256. - dstp[x] / 256.) / 4., 0.25)) *
			((srcp[x] / 256. - dstp[x] / 256.) / (_abs(srcp[x] / 256. - dstp[x] / 256.) + 1.001)) * 256.);
	_clamp(val, 16);
	dstp[x] = static_cast<uint16_t>(val);
}

static inline void _reverse_dif(const uint8_t *&srcp, uint8_t *&dstp, int x, double _amp, bool _linear) {
	long val = 0l;
	if (_linear)
		val = static_cast<long>((1. + _amp) * srcp[x] - dstp[x] * _amp);
	else
		val = static_cast<long>(srcp[x] +
			(_amp * 4. * pow(_abs(static_cast<double>(srcp[x]) - dstp[x]) / 4., 0.25)) *
			((static_cast<double>(srcp[x]) - dstp[x]) / (_abs(static_cast<double>(srcp[x]) - dstp[x]) + 1.001)));
	_clamp(val, 8);
	dstp[x] = static_cast<uint8_t>(val);
}

template<typename T1, typename T2>
static inline void _min_dif(const T1 *&srcp, T1 *&dstp, T1 *&dstp_a, T1 *&dstp_b, int x, double _amp) {
	T2 _dif_a = _abs(static_cast<T2>(dstp_a[x]) - srcp[x]);
	T2 _dif_b = _abs(static_cast<T2>(dstp_b[x]) - srcp[x]);
	dstp[x] = (_dif_a > _dif_b) ? dstp_b[x] : dstp_a[x];
}

template<typename T1, typename T2>
static void _apply_dif(const T1 *&srcp, int src_stride, T1 *&dstp, int dst_stride, int h, int w, T1 *&dstp_a, T1 *&dstp_b, double _amp, int mode, bool _linear) {
	const T1 *srctmp = srcp;
	T1 *dsttmp = dstp;
	T1 *dst_atmp = dstp_a;
	T1 *dst_btmp = dstp_b;
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			if (!mode)
				dstp[x] = srcp[x];
			else if (mode == 1) {
				dstp[x] = dstp_a[x];
				_reverse_dif(srcp, dstp, x, _amp, _linear);
			}
			else if (mode == 2) {
				dstp[x] = dstp_b[x];
				_reverse_dif(srcp, dstp, x, _amp, _linear);
			}
			else {
				_min_dif<T1, T2>(srcp, dstp, dstp_a, dstp_b, x, _amp);
				_reverse_dif(srcp, dstp, x, _amp, _linear);
			}
		}
		dstp += dst_stride;
		srcp += src_stride;
		dstp_a += dst_stride;
		dstp_b += dst_stride;
	}
	srcp = srctmp;
	dstp = dsttmp;
	dstp_a = dst_atmp;
	dstp_b = dst_btmp;
}

template<typename T1, typename T2>
static inline void _do_it_(const VSAPI *&vsapi, const VSFrameRef *&src, int src_stride, VSFrameRef *&dst, int dst_stride, VSFrameRef *&pad, int pad_stride, int h, int w, VSFrameRef *&dst_a, VSFrameRef *&dst_b, int plane, double _amp, int mode, bool _linear) {
	const T1 *srcp = reinterpret_cast<const T1 *>(vsapi->getReadPtr(src, plane));
	T1 *padp = reinterpret_cast<T1 *>(vsapi->getWritePtr(pad, plane));
	T1 *dstp = reinterpret_cast<T1 *>(vsapi->getWritePtr(dst, plane));
	T1 *dstp_a = reinterpret_cast<T1 *>(vsapi->getWritePtr(dst_a, plane));
	T1 *dstp_b = reinterpret_cast<T1 *>(vsapi->getWritePtr(dst_b, plane));
	if (!mode)
		_apply_dif<T1, T2>(srcp, src_stride, dstp, dst_stride, h, w, dstp_a, dstp_b, _amp, mode, _linear);
	else if (mode == 1) {
		_pad<T1, T2>(srcp, src_stride, padp, pad_stride, h, w);
		_process_plane_3x3(srcp, src_stride, dstp_a, dst_stride, h, w, true);
		_apply_dif<T1, T2>(srcp, src_stride, dstp, dst_stride, h, w, dstp_a, dstp_b, _amp, mode, _linear);
	}
	else if (mode == 2) {
		_pad<T1, T2>(srcp, src_stride, padp, pad_stride, h, w);
		_process_plane_3x3(srcp, src_stride, dstp_b, dst_stride, h, w, false);
		_apply_dif<T1, T2>(srcp, src_stride, dstp, dst_stride, h, w, dstp_a, dstp_b, _amp, mode, _linear);
	}
	else {
		_pad<T1, T2>(srcp, src_stride, padp, pad_stride, h, w);
		_process_plane_3x3(srcp, src_stride, dstp_a, dst_stride, h, w, true);
		_process_plane_3x3(srcp, src_stride, dstp_b, dst_stride, h, w, false);
		_apply_dif<T1, T2>(srcp, src_stride, dstp, dst_stride, h, w, dstp_a, dstp_b, _amp, mode, _linear);
	}
}

#endif