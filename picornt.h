#pragma once

int find_objects
		(
			float rs[], float cs[], float ss[], float qs[], int maxndetections,
			void* cascade, float angle,
			void* pixels, int nrows, int ncols, int ldim,
			float scalefactor, float stridefactor, float minsize, float maxsize
		);

int cluster_detections(float rs[], float cs[], float ss[], float qs[], int n);
