#include <cstdlib>
#include "Kernel.hpp"

struct MinSRPData {
	VSNodeRef *node;
	const VSVideoInfo *vi;
	double str[3];
	int mode[3];
};

static void VS_CC minsrpInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
	MinSRPData *d = reinterpret_cast<MinSRPData *>(*instanceData);
	vsapi->setVideoInfo(d->vi, 1, node);
}

static const VSFrameRef *VS_CC minsrpGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
	MinSRPData *d = reinterpret_cast<MinSRPData *>(*instanceData);
	d->vi = vsapi->getVideoInfo(d->node);
	if (activationReason == arInitial) {
		vsapi->requestFrameFilter(n, d->node, frameCtx);
	}
	else if (activationReason == arAllFramesReady) {
		const VSFrameRef *src = vsapi->getFrameFilter(n, d->node, frameCtx);
		const VSFormat *fi = d->vi->format;
		int height = vsapi->getFrameHeight(src, 0);
		int width = vsapi->getFrameWidth(src, 0);
		VSFrameRef *dst_a = vsapi->newVideoFrame(fi, width, height, src, core);
		VSFrameRef *dst_b = vsapi->newVideoFrame(fi, width, height, src, core);
		VSFrameRef *dst = vsapi->newVideoFrame(fi, width, height, src, core);
		VSFrameRef *pad = vsapi->newVideoFrame(fi, width + 2, height + 2, src, core);
		int plane;
		for (plane = 0; plane < fi->numPlanes; ++plane) {
			int src_stride = vsapi->getStride(src, plane) / d->vi->format->bytesPerSample;
			int dst_stride = vsapi->getStride(dst, plane) / d->vi->format->bytesPerSample;
			int pad_stride = vsapi->getStride(pad, plane) / d->vi->format->bytesPerSample;
			int h = vsapi->getFrameHeight(src, plane);
			int w = vsapi->getFrameWidth(src, plane);
			if (d->vi->format->bytesPerSample == 4)
				_do_it_<float, double>(vsapi, src, src_stride, dst, dst_stride, pad, pad_stride, h, w, dst_a, dst_b, plane, d->str[plane], d->mode[plane]);
			else if (d->vi->format->bytesPerSample > 1 && d->vi->format->bytesPerSample < 4)
				_do_it_<uint16_t, int32_t>(vsapi, src, src_stride, dst, dst_stride, pad, pad_stride, h, w, dst_a, dst_b, plane, d->str[plane], d->mode[plane]);
			else
				_do_it_<uint8_t, int32_t>(vsapi, src, src_stride, dst, dst_stride, pad, pad_stride, h, w, dst_a, dst_b, plane, d->str[plane], d->mode[plane]);
		}
		vsapi->freeFrame(src);
		vsapi->freeFrame(pad);
		vsapi->freeFrame(dst_a);
		vsapi->freeFrame(dst_b);
		return dst;
	}
	return nullptr;
}

static void VS_CC minsrpFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
	MinSRPData *d = reinterpret_cast<MinSRPData *>(instanceData);
	vsapi->freeNode(d->node);
	delete d;
}

static void VS_CC minsrpCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
	MinSRPData d;
	MinSRPData *data = nullptr;
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);
	const int m = vsapi->propNumElements(in, "mode");
	const int n = vsapi->propNumElements(in, "str");
	if (m > d.vi->format->numPlanes) {
		vsapi->setError(out, "MinSRP: number of modes specified must be equal to or fewer than the number of input planes");
		vsapi->freeNode(d.node);
		return;
	}
	if (n > d.vi->format->numPlanes) {
		vsapi->setError(out, "MinSRP: number of the specified elements in str array must be equal to or fewer than the number of input planes");
		vsapi->freeNode(d.node);
		return;
	}
	for (int i = 0; i < 3; ++i) {
		if (m == -1)
			d.mode[0] = d.mode[1] = d.mode[2] = 3;
		else
			if (i < m) {
				d.mode[i] = int64ToIntS(vsapi->propGetInt(in, "mode", i, nullptr));
				if (d.mode[i] < 0 || d.mode[i] > 3) {
					vsapi->setError(out, "MinSRP: invalid mode specified, only modes 0-3 supported");
					vsapi->freeNode(d.node);
					return;
				}
			}
			else
				d.mode[i] = d.mode[i - 1];
		if (n == -1)
			d.str[0] = d.str[1] = d.str[2] = 1.;
		else
			if (i < n)
				d.str[i] = vsapi->propGetFloat(in, "str", i, nullptr);
			else
				d.str[i] = d.str[i - 1];
	}
	if (!isConstantFormat(d.vi)) {
		vsapi->setError(out, "MinSRP: only input with constant format supported");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->subSamplingH || d.vi->format->subSamplingW) {
	        vsapi->setError(out, "MinSRP: 4:4:4 or gray input required!");
		vsapi->freeNode(d.node);
		return;
		}
	data = new MinSRPData;
	*data = d;
	vsapi->createFilter(in, out, "MinSRP", minsrpInit, minsrpGetFrame, minsrpFree, fmParallel, 0, data, core);
}

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin) {
	configFunc("com.sharpeners.minsrp", "minsrp", "VapourSynth MinSRP Filter", VAPOURSYNTH_API_VERSION, 1, plugin);
	registerFunc("Sharp",
		"clip:clip;"
		"str:float[]:opt;"
		"mode:int[]:opt;",
		minsrpCreate, 0, plugin);
}
