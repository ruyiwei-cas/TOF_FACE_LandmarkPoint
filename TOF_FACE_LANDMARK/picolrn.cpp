///*
// *	Copyright (c) 2013, Nenad Markus
// *	All rights reserved.
// *
// *	This is an implementation of the algorithm described in the following paper:
// *		N. Markus, M. Frljak, I. S. Pandzic, J. Ahlberg and R. Forchheimer,
// *		Object Detection with Pixel Intensity Comparisons Organized in Decision Trees,
// *		http://arxiv.org/abs/1305.4537
// *
// *	Redistribution and use of this program as source code or in binary form, with or without modifications, are permitted provided that the following conditions are met:
// *		1. Redistributions may not be sold, nor may they be used in a commercial product or activity without prior permission from the copyright holder (contact him at nenad.markus@fer.hr).
// *		2. Redistributions may not be used for military purposes.
// *		3. Any published work which utilizes this program shall include the reference to the paper available at http://arxiv.org/abs/1305.4537
// *		4. Redistributions must retain the above copyright notice and the reference to the algorithm on which the implementation is based on, this list of conditions and the following disclaimer.
// *
// *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// *	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// */
//
//#include <stdio.h>
//#include <malloc.h>
//#include <math.h>
//#include <stdint.h>
//
//// hyperparameters
//#define NRANDS 1024
//
///*
//	auxiliary stuff
//*/
//
//#define MAX(a, b) ((a)>(b)?(a):(b))
//#define MIN(a, b) ((a)<(b)?(a):(b))
//#define SQR(x) ((x)*(x))
//
////这段代码值得借鉴
///*
//	portable time function
//*/
//
//#ifdef __GNUC__
//#include <time.h>
//float getticks()
//{
//	struct timespec ts;
//
//	if(clock_gettime(CLOCK_MONOTONIC, &ts) < 0)
//		return -1.0f;
//
//	return ts.tv_sec + 1e-9f*ts.tv_nsec;
//}
//#else
//#include <windows.h>
//float getticks()
//{
//	static double freq = -1.0;
//	LARGE_INTEGER lint;
//
//	if(freq < 0.0)
//	{
//		if(!QueryPerformanceFrequency(&lint))
//			return -1.0f;
//
//		freq = lint.QuadPart;
//	}
//
//	if(!QueryPerformanceCounter(&lint))
//		return -1.0f;
//
//	return (float)( lint.QuadPart/freq );
//}
//#endif
//
///*
//	随机树生成器
//	multiply with carry PRNG
//*/
//uint32_t mwcrand_r(uint64_t* state)
//{
//	uint32_t* m;
//
//	//
//	m = (uint32_t*)state;
//
//	// bad state?
//	if(m[0] == 0)
//		m[0] = 0xAAAA;
//
//	if(m[1] == 0)
//		m[1] = 0xBBBB;
//
//	// mutate state
//	m[0] = 36969 * (m[0] & 65535) + (m[0] >> 16);
//	m[1] = 18000 * (m[1] & 65535) + (m[1] >> 16);
//
//	// output
//	return (m[0] << 16) + m[1];
//}
//
//uint64_t prngglobal = 0x12345678000fffffLL;
//
//void smwcrand(uint32_t seed)
//{
//	prngglobal = 0x12345678000fffffLL*seed;
//}
//
//uint32_t mwcrand()
//{
//	return mwcrand_r(&prngglobal);
//}
//
///*
//	
//*/
//
//#define MAX_N 2000000
//
//int N = 0;
////存储图像，最大2000000个
//uint8_t* ppixels[MAX_N];
//int pdims[MAX_N][2]; // (nrows, ncols) 每张图的长宽
//
//int nbackground = 0;
//// 0,2,4，5，6 表示在ppixels中第0，2，4，5，6为负样本
//int background[MAX_N]; // i
////样本的总共数目
//int nobjects = 0;
//// 1,3,9 表示在ppixels中第1，3，9为正样本
//int objects[MAX_N][4]; // (r, c, s, i) rcs row col size的意思，每个样本都有这4个数据，在这个数据集合中，size的大小是一个正方形，col row是指top，left的位置
//
//int load_image(uint8_t* pixels[], int* nrows, int* ncols, FILE* file)
//{
//	/*
//	- loads an 8-bit grey image saved in the <RID> file format
//	- <RID> file contents:
//		- a 32-bit signed integer h (image height)
//		- a 32-bit signed integer w (image width)
//		- an array of w*h unsigned bytes representing pixel intensities
//	*/
//
//	//
//	if(fread(nrows, sizeof(int), 1, file) != 1)
//		return 0;
//
//	if(fread(ncols, sizeof(int), 1, file) != 1)
//		return 0;
//
//	//
//	*pixels = (uint8_t*)malloc(*nrows**ncols*sizeof(uint8_t));
//
//	if(!*pixels)
//		return 0;
//
//	// read pixels
//	if(fread(*pixels, sizeof(uint8_t), *nrows**ncols, file) != *nrows**ncols)
//		return 0;
//
//	// we're done
//	return 1;
//}
//
//int load_training_data(char* path)
//{
//	FILE* file;
//
//	//
//	file = fopen(path, "rb");
//
//	if(!file)
//		return 0;
//	//
//	N = 0;
//
//	nbackground = 0;
//	nobjects = 0;
//
//	while( load_image(&ppixels[N], &pdims[N][0], &pdims[N][1], file) )
//	{
//		int i, n;
//
//		//文件最后是一个label
//		if(fread(&n, sizeof(int), 1, file) != 1)
//			return 1;
//
//		//负样本是0
//		if(!n)
//		{
//			//标记background第n个对应总体数据的第N个
//			background[nbackground] = N;
//			++nbackground;
//		}
//		else
//		{
//			//正样本
//			for(i=0; i<n; ++i)
//			{
//				//如果是正样本，那么它还有 r c s需要读取
//				fread(&objects[nobjects][0], sizeof(int), 1, file); // r 
//				fread(&objects[nobjects][1], sizeof(int), 1, file); // c
//				fread(&objects[nobjects][2], sizeof(int), 1, file); // s
//
//				objects[nobjects][3] = N; // i 标记它是总体样本中的第几个样本
//				//
//				++nobjects;
//			}
//		}
//		//
//		++N;
//	}
//
//	//
//	return 1;
//}
//
///*
//	regression trees
//*/
//
//int bintest(int32_t tcode, int r, int c, int s, int iind)
//{
//	//
//	int r1, c1, r2, c2;
//	int8_t* p = (int8_t*)&tcode;
//
//	//
//	r1 = (256*r + p[0]*s)/256;
//	c1 = (256*c + p[1]*s)/256;
//
//	r2 = (256*r + p[2]*s)/256;
//	c2 = (256*c + p[3]*s)/256;
//
//	//
//	r1 = MIN(MAX(0, r1), pdims[iind][0]-1);
//	c1 = MIN(MAX(0, c1), pdims[iind][1]-1);
//
//	r2 = MIN(MAX(0, r2), pdims[iind][0]-1);
//	c2 = MIN(MAX(0, c2), pdims[iind][1]-1);
//
//	//
//	return ppixels[iind][r1*pdims[iind][1]+c1]<=ppixels[iind][r2*pdims[iind][1]+c2];
//}
//
//float get_split_error(int32_t tcode, float tvals[], int rs[], int cs[], int ss[], int iinds[], double ws[], int inds[], int indsnum)
//{
//	int i, j;
//
//	double wsum, wsum0, wsum1;
//	double wtvalsum0, wtvalsumsqr0, wtvalsum1, wtvalsumsqr1;
//
//	double wmse0, wmse1;
//
//	//
//	wsum = wsum0 = wsum1 = wtvalsum0 = wtvalsum1 = wtvalsumsqr0 = wtvalsumsqr1 = 0.0;
//
//	for(i=0; i<indsnum; ++i)
//	{
//		if( bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], iinds[inds[i]]) )
//		{
//			wsum1 += ws[inds[i]];
//			wtvalsum1 += ws[inds[i]]*tvals[inds[i]];
//			wtvalsumsqr1 += ws[inds[i]]*SQR(tvals[inds[i]]);
//		}
//		else
//		{
//			wsum0 += ws[inds[i]];
//			wtvalsum0 += ws[inds[i]]*tvals[inds[i]];
//			wtvalsumsqr0 += ws[inds[i]]*SQR(tvals[inds[i]]);
//		}
//
//		wsum += ws[inds[i]];//总共的weight
//	}
//
//	//
//	wmse0 = wtvalsumsqr0 - SQR(wtvalsum0)/wsum0;
//	wmse1 = wtvalsumsqr1 - SQR(wtvalsum1)/wsum1;
//
//	//
//	return (float)( (wmse0 + wmse1)/wsum );
//}
//
//int split_training_data(int32_t tcode, float tvals[], int rs[], int cs[], int ss[], int iinds[], double ws[], int inds[], int ninds)
//{
//	int stop;
//	int i, j;
//
//	int n0;
//
//	//
//	stop = 0;
//
//	i = 0;
//	j = ninds - 1;
//
//	while(!stop)
//	{
//		//交换正负样本的位置
//		while( !bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], iinds[inds[i]]) )
//		{
//			if( i==j )
//				break;
//			else
//				++i;
//		}
//
//		while( bintest(tcode, rs[inds[j]], cs[inds[j]], ss[inds[j]], iinds[inds[j]]) )
//		{
//			if( i==j )
//				break;
//			else
//				--j;
//		}
//
//		//
//		if( i==j )
//			stop = 1;
//		else
//		{
//			// swap
//			inds[i] = inds[i] ^ inds[j];
//			inds[j] = inds[i] ^ inds[j];
//			inds[i] = inds[i] ^ inds[j];
//		}
//	}
//
//	//
//	n0 = 0;
//
//	for (i = 0; i < ninds; ++i){
//		if (!bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], iinds[inds[i]])){
//			++n0;
//		}
//	}
//	//
//	return n0;
//}
//
//int32_t get_random_tcode(int8_t* bbox)
//{
//	int32_t tcode;
//	int8_t* p;
//
//	//
//	p = (int8_t*)&tcode;
//
//	//
//	p[0] = bbox[0] + mwcrand()%(bbox[1]-bbox[0]+1);
//	p[1] = bbox[2] + mwcrand()%(bbox[3]-bbox[2]+1);
//
//	p[2] = bbox[0] + mwcrand()%(bbox[1]-bbox[0]+1);
//	p[3] = bbox[2] + mwcrand()%(bbox[3]-bbox[2]+1);
//
//	//
//	return tcode;
//}
////d是表示树当前的深度，maxd表示树的最大深度
//int grow_subtree(int32_t tcodes[], float lut[], int nodeidx, int d, int maxd, float tvals[], int rs[], int cs[], int ss[], int iinds[], double ws[], int inds[], int ninds, int8_t* bbox)
//{
//	int i, nrands;
//
//	int32_t tmptcodes[2048];
//	float es[2048], e;
//
//	int n0;
//
//	//到达最大深度
//	if(d == maxd)
//	{
//		int lutidx;
//		double tvalaccum, wsum;
//
//		//
//		lutidx = nodeidx - ((1<<maxd)-1);
//
//		// compute output: a simple average
//		tvalaccum = 0.0;
//		wsum = 0.0;
//
//		for(i=0; i<ninds; ++i)
//		{
//			tvalaccum += ws[inds[i]]*tvals[inds[i]];
//			wsum += ws[inds[i]];
//		}
//
//		if(wsum == 0.0)
//			lut[lutidx] = 0.0f;
//		else
//			lut[lutidx] = (float)( tvalaccum/wsum );
//
//		//
//		return 1;
//	}
//	else if(ninds <= 1) //样本数据只有1个
//	{
//		//
//		tcodes[nodeidx] = 0;
//
//		//
//		grow_subtree(tcodes, lut, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, ss, iinds, ws, inds, ninds, bbox);
//		grow_subtree(tcodes, lut, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, ss, iinds, ws, inds, ninds, bbox);
//
//		return 1;
//	}
//
//	// generate binary test codes
//	//测试的测试是
//	nrands = NRANDS;
//
//	//随机获取4个点
//	for (i = 0; i < nrands; ++i){
//		tmptcodes[i] = get_random_tcode(bbox);
//	}
//	//
//	#pragma omp parallel for
//	for (i = 0; i < nrands; ++i){
//		es[i] = get_split_error(tmptcodes[i], tvals, rs, cs, ss, iinds, ws, inds, ninds);
//	}
//	//
//	e = es[0];
//	tcodes[nodeidx] = tmptcodes[0];
//
//	for (i = 1; i<nrands; ++i){
//		if (e > es[i])
//		{
//			e = es[i];
//			tcodes[nodeidx] = tmptcodes[i];
//		}
//	}
//	//
//	n0 = split_training_data(tcodes[nodeidx], tvals, rs, cs, ss, iinds, ws, inds, ninds);
//
//	//把整个样本分成两个部分，
//	grow_subtree(tcodes, lut, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, ss, iinds, ws, &inds[0], n0, bbox);
//	grow_subtree(tcodes, lut, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, ss, iinds, ws, &inds[n0], ninds-n0, bbox);
//
//	//
//	return 1;
//}
//
//int grow_rtree(int32_t tcodes[], float lut[], int d, float tvals[], int rs[], int cs[], int ss[], int iinds[], double ws[], int n, int8_t* bbox)
//{
//	int i;
//	int* inds;
//
//	//n是指样本的个数，
//	inds = (int*)malloc(n*sizeof(int));
//	//每个样本的编号
//	for(i=0; i<n; ++i)
//		inds[i] = i;
//
//	//d表示数的深度
//	if(!grow_subtree(tcodes, lut, 0, 0, d, tvals, rs, cs, ss, iinds, ws, inds, n, bbox))
//	{
//		free(inds);
//		return 0;
//	}
//	else
//	{
//		free(inds);
//		return 1;
//	}
//}
//
///*
//	
//*/
//
//int32_t version = 3;
//
//int tdepth;
//int ntrees=0;
//
//int8_t bbox[4]; // (r_min, r_max, c_min, c_max)
//
////整个模型中，最多是4096个树，并且每个树的深度最大是9层
//int32_t tcodes[4096][1024];
////luts是指叶子节点上的值
//float luts[4096][1024];
//
//float thresholds[4096];
//
///*
//	
//*/
//
//int load_cascade_from_file(const char* path)
//{
//	int i;
//	FILE* file;
//
//	//
//	file = fopen(path, "rb");
//
//	if(!file)
//		return 0;
//
//	//
//	fread(&version, sizeof(int32_t), 1, file);
//	fread(&bbox[0], sizeof(int8_t), 4, file);
//	fread(&tdepth, sizeof(int), 1, file);
//	fread(&ntrees, sizeof(int), 1, file);
//
//	for(i=0; i<ntrees; ++i)
//	{
//		//
//		fread(&tcodes[i][0], sizeof(int32_t), (1<<tdepth)-1, file);
//		fread(&luts[i][0], sizeof(float), 1<<tdepth, file);
//		fread(&thresholds[i], sizeof(float), 1, file);
//	}
//
//	//
//	fclose(file);
//
//	//
//	return 1;
//}
//
//int save_cascade_to_file(const char* path)
//{
//	int i;
//	FILE* file;
//
//	//
//	file = fopen(path, "wb");
//
//	if(!file)
//		return 0;
//
//	//
//	fwrite(&version, sizeof(int32_t), 1, file);
//	fwrite(&bbox[0], sizeof(int8_t), 4, file);
//	fwrite(&tdepth, sizeof(int), 1, file);
//	fwrite(&ntrees, sizeof(int), 1, file);
//
//	for(i=0; i<ntrees; ++i)
//	{
//		//
//		fwrite(&tcodes[i][0], sizeof(int32_t), (1<<tdepth)-1, file);
//		fwrite(&luts[i][0], sizeof(float), 1<<tdepth, file);
//		fwrite(&thresholds[i], sizeof(float), 1, file);
//	}
//
//	//
//	fclose(file);
//
//	//
//	return 1;
//}
//
///*
//	
//*/
//
//float get_tree_output(int i, int r, int c, int s, int iind)
//{
//	int idx, j;
//
//	//
//	idx = 1;
//
//	for(j=0; j<tdepth; ++j)
//		idx = 2*idx + bintest(tcodes[i][idx-1], r, c, s, iind);
//
//	//
//	return luts[i][idx - (1<<tdepth)];
//}
//
//int classify_region(float* o, int r, int c, int s, int iind)
//{
//	int i, sr, sc;
//
//	//默认分为正样本，当正样本导入的时候返回正，当负样本导入的时候返回负，外面可以判断
//	if(!ntrees)
//		return 1;
//
//	//
//	*o = 0.0f;
//
//	for(i=0; i<ntrees; ++i)
//	{
//		//
//		*o += get_tree_output(i, r, c, s, iind);
//
//		//
//		if(*o <= thresholds[i])
//			return -1;
//	}
//
//	//
//	return 1;
//}
//
///*
//np：正样本的数目
//nn :负样本的数目
//mintpr：true positve ratio
//maxfp: false positve ratio 
//maxn trees:树的数目
//*/
//int learn_new_stage(float mintpr, float maxfpr, int maxntrees, float tvals[], int rs[], int cs[], int ss[], int iinds[], float os[], int np, int nn)
//{
//	int i;
//
//	double* ws;
//	double wsum;
//
//	float threshold, tpr, fpr;
//
//	//
//	printf("* learning a new stage ...\n");
//
//	//每个样本的权值
//	ws = (double*)malloc((np+nn)*sizeof(double));
//
//	//最大样本数
//	maxntrees = ntrees + maxntrees;
//	fpr = 1.0f;//被错误的判读为真的样本
//
//	//两个条件，一个是树的数目是小于max， fpr大于最大的fp
//	while(ntrees<maxntrees && fpr>maxfpr)
//	{
//		float t;
//		int numtps, numfps;
//
//		//
//		t = getticks();
//
//		// compute weights ...
//		wsum = 0.0;
//
//		for(i=0; i<np+nn; ++i)
//		{
//			if(tvals[i] > 0)
//				ws[i] = exp(-1.0*os[i])/np;
//			else
//				ws[i] = exp(+1.0*os[i])/nn;
//
//			wsum += ws[i]; //@COMPARE_ERROR 我的做法是直接对各个进行求解均值
//		}
//
//		for(i=0; i<np+nn; ++i) //更新每个样本的
//			ws[i] /= wsum;
//
//		// grow a tree ...
//		grow_rtree(tcodes[ntrees], luts[ntrees], tdepth, tvals, rs, cs, ss, iinds, ws, np+nn, bbox);
//
//		thresholds[ntrees] = -1337.0f;
//
//		++ntrees;
//
//		// update outputs ...
//		for(i=0; i<np+nn; ++i)
//		{
//			float o;
//
//			//
//			o = get_tree_output(ntrees-1, rs[i], cs[i], ss[i], iinds[i]);
//
//			//
//			os[i] += o;
//		}
//
//		// get threshold ... 每个tree都存在一个threshold，这个threshold是用于做为cascade的效果
//		threshold = 5.0f;
//		do
//		{
//			//
//			threshold -= 0.005f;
//
//			numtps = 0;
//			numfps = 0;
//
//			//
//			for(i=0; i<np+nn; ++i)
//			{
//				if( tvals[i]>0 && os[i]>threshold)
//					++numtps;
//				if(	tvals[i]<0 && os[i]>threshold)
//					++numfps;
//			}
//
//			//
//			tpr = numtps/(float)np;
//			fpr = numfps/(float)nn;
//		}
//		while(tpr<mintpr);
//
//		printf("	** tree %d (%d [s]) ... stage tpr=%f, stage fpr=%f\n", ntrees, (int)(getticks()-t), tpr, fpr);
//		fflush(stdout);
//	}
//
//	//最后一个threshold 才有作用，之前的threshold应该没有作用
//	thresholds[ntrees-1] = threshold;
//
//	printf("	** threshold set to %f\n", threshold);
//
//	//
//	free(ws);
//
//	//
//	return 1;
//}
//
///*tvals=label rs=row cs=col ss=size iinds=index np=positive numbe nn=negtive number*/
//float sample_training_data(float tvals[], int rs[], int cs[], int ss[], int iinds[], float os[], int* np, int* nn)
//{
//	int i, n;
//
//	int64_t nw;
//	float etpr, efpr;
//
//	int t;
//
//	#define NUMPRNGS 1024
//	static int prngsinitialized = 0;
//	static uint64_t prngs[NUMPRNGS];
//
//	int stop;
//
//	//
//	t = getticks();
//
//	//
//	n = 0;
//
//	/*
//		object samples
//	*/
//	//获取所有正确分类的正样本
//	for (i = 0; i < nobjects; ++i){
//		if( classify_region(&os[n], objects[i][0], objects[i][1], objects[i][2], objects[i][3]) == 1 )
//		{
//			//
//			rs[n] = objects[i][0];
//			cs[n] = objects[i][1];
//			ss[n] = objects[i][2];
//
//
//			iinds[n] = objects[i][3];
//
//			tvals[n] = +1;
//
//			//
//			++n;
//		}
//	}
//	*np = n;
//
//	/*
//		non-object samples
//	*/
//
//	if(!prngsinitialized)
//	{
//		// initialize a PRNG for each thread
//		for (i = 0; i < NUMPRNGS; ++i){
//			prngs[i] = 0xFFFF*mwcrand() + 0xFFFF1234FFFF0001LL*mwcrand();
//		}
//		prngsinitialized = 1;
//	}
//
//	//
//	nw = 0;
//	*nn = 0;
//
//	stop = 0;
//
//	if(nbackground)
//	{
//		#pragma omp parallel
//		{
//			int thid;
//
//			//
//			thid = omp_get_thread_num();
//
//			while(!stop)
//			{
//				/*
//					data mine hard negatives
//				*/
//
//				float o;
//				int iind, s, r, c, nrows, ncols;
//				uint8_t* pixels;
//
//				//获取background的样本
//				iind = background[ mwcrand_r(&prngs[thid])%nbackground ];
//
//				//随机在object的图像上获取一个位置
//				r = mwcrand_r(&prngs[thid])%pdims[iind][0];
//				c = mwcrand_r(&prngs[thid])%pdims[iind][1];
//				//随机获取一个size
//				s = objects[mwcrand_r(&prngs[thid])%nobjects][2]; // sample the size of a random object in the pool
//
//				//在背景数据上随机获取一个正方形做为一个负样本，当这个负样本检测我正样本的时候，那么，这是一个false postive数据
//				if( classify_region(&o, r, c, s, iind) == 1 )
//				{
//					//we have a false positive ...//把一个负样本认为是一个正样本
//					#pragma omp critical
//					{
//						if(*nn<*np)//获取一定数量的false postive 样本，
//						{
//							rs[n] = r;
//							cs[n] = c;
//							ss[n] = s;
//
//							iinds[n] = iind;
//
//							os[n] = o;
//
//							tvals[n] = -1;
//
//							//
//							++n;
//							++*nn;
//						}
//						else{
//							stop = 1; 
//						}
//					}
//				}
//
//				if(!stop)
//				{
//					#pragma omp atomic
//					++nw;
//				}
//			}
//		}
//	}
//	else{
//		nw = 1;
//	}
//	/*
//		print the estimated true positive and false positive rates
//	*/
//
//	etpr = *np/(float)nobjects;
//	efpr = (float)( *nn/(double)nw );
//
//	printf("* sampling finished ...\n");
//	printf("	** elapsed time: %d\n", (int)(getticks()-t));
//	printf("	** cascade TPR=%.8f\n", etpr);
//	printf("	** cascade FPR=%.8f (%d/%lld)\n", efpr, *nn, (long long int)nw);
//
//	/*
//		
//	*/
//
//	return efpr;
//}
//
///*
//	
//*/
//
//static int rs[2*MAX_N];
//static int cs[2*MAX_N];
//static int ss[2*MAX_N];
//static int iinds[2*MAX_N];
//static float tvals[2*MAX_N];
//static float os[2*MAX_N];
//
//int learn_with_default_parameters(char* trdata, char* dst)
//{
//	int i, np, nn;
//	float fpr;
//
//	//
//	if(!load_training_data(trdata))
//	{
//		printf("* cannot load training data ...\n");
//		return 0;
//	}
//
//	//
//	bbox[0] = -127;
//	bbox[1] = +127;
//	bbox[2] = -127;
//	bbox[3] = +127;
//	//随机树的层数
//	tdepth = 5;
//
//	if(!save_cascade_to_file(dst))
//			return 0;
//
//	//
//	sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
//	learn_new_stage(0.9800f, 0.5f, 4, tvals, rs, cs, ss, iinds, os, np, nn);
//	save_cascade_to_file(dst);
//
//	printf("\n");
//
//	sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
//	learn_new_stage(0.9850f, 0.5f, 8, tvals, rs, cs, ss, iinds, os, np, nn);
//	save_cascade_to_file(dst);
//
//	printf("\n");
//
//	sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
//	learn_new_stage(0.9900f, 0.5f, 16, tvals, rs, cs, ss, iinds, os, np, nn);
//	save_cascade_to_file(dst);
//
//	printf("\n");
//
//	sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
//	learn_new_stage(0.9950f, 0.5f, 32, tvals, rs, cs, ss, iinds, os, np, nn);
//	save_cascade_to_file(dst);
//
//	printf("\n");
//
//	//
//	while(sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn) > 1e-6f)
//	{
//		learn_new_stage(0.9975f, 0.5f, 64, tvals, rs, cs, ss, iinds, os, np, nn);
//		save_cascade_to_file(dst);
//
//		printf("\n");
//	}
//
//	//
//	printf("* target FPR achieved ... terminating the learning process ...\n");
//}
//
///*
//	
//*/
//
//const char* howto()
//{
//	return
//		"TODO\n"
//	;
//}
//
//int main(int argc, char* argv[])
//{
//	// initialize the PRNG
//	smwcrand(time(0));
//
//	//
//	if(argc == 3)
//	{
//		learn_with_default_parameters(argv[1], argv[2]);
//	}
//	else if(argc == 7)
//	{
//		int dummy;
//
//		//
//		sscanf(argv[1], "%d", &dummy); bbox[0] = dummy;
//		sscanf(argv[2], "%d", &dummy); bbox[1] = dummy;
//		sscanf(argv[3], "%d", &dummy); bbox[2] = dummy;
//		sscanf(argv[4], "%d", &dummy); bbox[3] = dummy;
//
//		//
//		sscanf(argv[5], "%d", &tdepth);
//
//		//
//		ntrees = 0;
//
//		//
//		if(!save_cascade_to_file(argv[6]))
//			return 0;
//
//		//
//		printf("* initializing:\n");
//		printf("	** bbox = (%d, %d, %d, %d)\n", bbox[0], bbox[1], bbox[2], bbox[3]);
//		printf("	** tdepth = %d\n", tdepth);
//
//		//
//		return 0;
//	}
//	else if(argc == 7)
//	{
//		float tpr, fpr;
//		int ntrees, np, nn;
//
//		//
//		if(!load_cascade_from_file(argv[1]))
//		{
//			printf("* cannot load a cascade from '%s'\n", argv[1]);
//			return 1;
//		}
//
//		if(!load_training_data(argv[2]))
//		{
//			printf("* cannot load the training data from '%s'\n", argv[2]);
//			return 1;
//		}
//
//		//
//		sscanf(argv[3], "%f", &tpr);
//		sscanf(argv[4], "%f", &fpr);
//		sscanf(argv[5], "%d", &ntrees);
//
//		//
//		sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
//		learn_new_stage(tpr, fpr, ntrees, tvals, rs, cs, ss, iinds, os, np, nn);
//
//		//
//		if(!save_cascade_to_file(argv[6]))
//			return 1;
//	}
//	else
//	{
//		printf("%s", howto());
//		return 0;
//	}
//
//	//
//	return 0;
//}
