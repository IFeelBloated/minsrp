#ifndef KERNEL_HPP
#define KERNEL_HPP

#define _JUMP \
        dstp += dst_stride; \
        srcp += src_stride
#define _PERFORM_GAUSS \
        _process_plane_3x3(srcp, src_stride, dstp_a, dst_stride, h, w, true)
#define _PERFORM_MEDIAN \
        _process_plane_3x3(srcp, src_stride, dstp_b, dst_stride, h, w, false)
#define _AMPLIFY_DIF \
        _reverse_dif(srcp, dstp, x, _amp, _linear)
#define _PERFORM_PAD \
        _pad<T1, T2>(srcp, src_stride, padp, pad_stride, h, w)

#define _SECURE_SCL 0

#include <algorithm>
#include "VapourSynth.h"
#include "VSHelper.h"
using namespace std;

static volatile uint64_t _ptr_sucker = 0xFFFFFFFFFFFFFFFF;
static long val = 0;

template<typename T>
static inline T _abs(T val) {
	return (val > 0) ? val : -val;
}

static inline void _clamp(int bits) {
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

static void _gauss_3x3(float *dstp, int x) {
	static float * const a = reinterpret_cast<float *>(_ptr_sucker);
	dstp[x] = static_cast<float>(
		(4. * a[0] +
			2. * (static_cast<double>(a[2]) + a[4] + a[5] + a[7]) +
			a[1] + a[3] + a[6] + a[8]) / 16.);
}

static void _gauss_3x3(uint16_t *dstp, int x) {
	static uint16_t * const a = reinterpret_cast<uint16_t *>(_ptr_sucker);
	val = ((a[0] << 2l) +
		((static_cast<uint32_t>(a[2]) + a[4] + a[5] + a[7]) << 1l) +
		a[1] + a[3] + a[6] + a[8] + 8l) >> 4l;
	_clamp(16);
	dstp[x] = static_cast<uint16_t>(val);
}

static void _gauss_3x3(uint8_t *dstp, int x) {
	static uint8_t * const a = reinterpret_cast<uint8_t *>(_ptr_sucker);
	val = ((a[0] << 2l) +
		((static_cast<uint32_t>(a[2]) + a[4] + a[5] + a[7]) << 1l) +
		a[1] + a[3] + a[6] + a[8] + 8l) >> 4l;
	_clamp(8);
	dstp[x] = static_cast<uint8_t>(val);
}

template<typename T>
static void _median_3x3(T *dstp, int x) {
	static T * const a = reinterpret_cast<T *>(_ptr_sucker);
	nth_element(a, a + 4, a + 9);
	dstp[x] = a[4];
}

template<typename T>
static void _process_plane_3x3(const T *srcp, int src_stride, T *dstp, int dst_stride, int h, int w, bool gauss) {
	void(*_func)(T *, int) = nullptr;
	_declspec(align(16)) static T a[9];
	_ptr_sucker = reinterpret_cast<uint64_t>(a);
	if (gauss)
		_func = _gauss_3x3;
	else
		_func = _median_3x3<T>;
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			a[0] = *(srcp + x);
			a[1] = *(srcp + x - src_stride - 1);
			a[2] = *(srcp + x - src_stride);
			a[3] = *(srcp + x - src_stride + 1);
			a[4] = *(srcp + x - 1);
			a[5] = *(srcp + x + 1);
			a[6] = *(srcp + x + src_stride - 1);
			a[7] = *(srcp + x + src_stride);
			a[8] = *(srcp + x + src_stride + 1);
			_func(dstp, x);
		}
		_JUMP;
	}
}

static void _reverse_dif(const float *srcp, float *dstp, int x, double _amp, bool _linear) {
	if (_linear)
		dstp[x] = static_cast<float>((1. + _amp) * srcp[x] - dstp[x] * _amp);
	else
		dstp[x] = static_cast<float>(srcp[x] +
			(_amp * 4. * pow(_abs(srcp[x] * 256. - dstp[x] * 256.) / 4., 0.25)) *
			((srcp[x] * 256. - dstp[x] * 256.) / (_abs(srcp[x] * 256. - dstp[x] * 256.) + 1.001)) / 256.);
}

static void _reverse_dif(const uint16_t *srcp, uint16_t *dstp, int x, double _amp, bool _linear) {
	if (_linear)
		val = static_cast<long>((1. + _amp) * srcp[x] - dstp[x] * _amp);
	else
		val = static_cast<long>(srcp[x] +
			(_amp * 4. * pow(_abs(srcp[x] / 256. - dstp[x] / 256.) / 4., 0.25)) *
			((srcp[x] / 256. - dstp[x] / 256.) / (_abs(srcp[x] / 256. - dstp[x] / 256.) + 1.001)) * 256.);
	_clamp(16);
	dstp[x] = static_cast<uint16_t>(val);
}

static void _reverse_dif(const uint8_t *srcp, uint8_t *dstp, int x, double _amp, bool _linear) {
	if (_linear)
		val = static_cast<long>((1. + _amp) * srcp[x] - dstp[x] * _amp);
	else
		val = static_cast<long>(srcp[x] +
			(_amp * 4. * pow(_abs(static_cast<double>(srcp[x]) - dstp[x]) / 4., 0.25)) *
			((static_cast<double>(srcp[x]) - dstp[x]) / (_abs(static_cast<double>(srcp[x]) - dstp[x]) + 1.001)));
	_clamp(8);
	dstp[x] = static_cast<uint8_t>(val);
}

template<typename T1, typename T2>
static void _do_it_(const VSAPI *vsapi, const VSFrameRef *src, int src_stride, VSFrameRef *dst, int dst_stride, VSFrameRef *pad, int pad_stride, int h, int w, VSFrameRef *dst_a, VSFrameRef *dst_b, int plane, double _amp, int mode, bool _linear) {
	const T1 *srcp = reinterpret_cast<const T1 *>(vsapi->getReadPtr(src, plane));
	T1 *padp = reinterpret_cast<T1 *>(vsapi->getWritePtr(pad, plane));
	T1 *dstp = reinterpret_cast<T1 *>(vsapi->getWritePtr(dst, plane));
	T1 *dstp_a = reinterpret_cast<T1 *>(vsapi->getWritePtr(dst_a, plane));
	T1 *dstp_b = reinterpret_cast<T1 *>(vsapi->getWritePtr(dst_b, plane));
	if (!mode)
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x)
				dstp[x] = srcp[x];
			_JUMP;
		}
	else if (mode == 1) {
		_PERFORM_PAD;
		_PERFORM_GAUSS;
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				dstp[x] = dstp_a[x];
				_AMPLIFY_DIF;
			}
			_JUMP;
			dstp_a += dst_stride;
		}
	}
	else if (mode == 2) {
		_PERFORM_PAD;
		_PERFORM_MEDIAN;
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				dstp[x] = dstp_b[x];
				_AMPLIFY_DIF;
			}
			_JUMP;
			dstp_b += dst_stride;
		}
	}
	else {
		_PERFORM_PAD;
		_PERFORM_GAUSS;
		_PERFORM_MEDIAN;
		T2 _dif_a = 0;
		T2 _dif_b = 0;
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				_dif_a = _abs(static_cast<T2>(dstp_a[x]) - srcp[x]);
				_dif_b = _abs(static_cast<T2>(dstp_b[x]) - srcp[x]);
				dstp[x] = (_dif_a > _dif_b) ? dstp_b[x] : dstp_a[x];
				_AMPLIFY_DIF;
			}
			_JUMP;
			dstp_a += dst_stride;
			dstp_b += dst_stride;
		}
	}
}

#endif
