#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double dist(int ay, int ax, int by, int bx) {
	return sqrt(pow(ax - bx, 2.0) + pow(ay - by, 2.0));
}

int max(int arr[], int len) {
	int max = arr[0];
	for (int i = 1; i < len; i++)
		if (arr[i] > max)
			max = arr[i];
	return max;
}

int min(int arr[], int len) {
	int min = arr[0];
	for (int i = 1; i < len; i++)
		if (arr[i] < min)
			min = arr[i];
	return min;
}

int* threelargest(int accumulator[], int acc_size, double radius, int xs[], int ys[]) {	 //TEMPORARY -- prints 3 most promising annuli for a given accumulator matrix
	int a, b, c, ai, bi, ci; //see int *hough(...) for more on this
	a = b = c = 0;
	for (int i = 0; i < acc_size; i++) {
		if (accumulator[i] > accumulator[a]) {
			c = b;
			b = a;
			a = i;
		}
		else if (accumulator[i] > accumulator[b]) {
			c = b;
			b = i;
		}
		else if (accumulator[i] > accumulator[c])
			c = i;
	}
	printf("(%d, %d), r=%f: %d\n", xs[a], ys[a], radius, accumulator[a]);
	printf("(%d, %d), r=%f: %d\n", xs[b], ys[b], radius, accumulator[b]);
	printf("(%d, %d), r=%f: %d\n", xs[c], ys[c], radius, accumulator[c]);

	int *ret = malloc(3 * sizeof(int));
	ret[0] = accumulator[a];
	ret[1] = xs[a];
	ret[2] = ys[a];
	return ret;
}

int *hough(int xs[], int ys[], int photonCount) {
	int rows = max(ys, photonCount) - min(ys, photonCount);
	int cols = max(xs, photonCount) - min(xs, photonCount);
	int minRadius = (rows < cols ? rows : cols) / 4.0; //the smallest radius size you want to test
	int maxRadius = (rows < cols ? rows : cols) / 1.5;
	int nRadii = 20;	//number of radii you want to test
	double fuzzy = 0.3;	//percent of radius leeway you want to give the photons

	int maxAcc = 0, argmaxX = 0, argmaxY = 0, argmaxRad = 0;

	for (double radius = minRadius; radius <= maxRadius; radius += (maxRadius - minRadius) / (double)nRadii) {
		printf("Generating Accumulator for r=%f\n", radius); //allocates accumulator matrix
		int *accumulator = malloc(photonCount * sizeof(int));

		for (int centralPhoton = 0; centralPhoton < photonCount; centralPhoton++) { //iterates over photons in the image to check distances of other photons
			int centerX = xs[centralPhoton];
			int centerY = ys[centralPhoton];
			int acc = 0;
			for (int testPhoton = 0; testPhoton < photonCount; testPhoton++) {
				int photonX = xs[testPhoton];
				int photonY = ys[testPhoton];
				if (centralPhoton != testPhoton) { //checks if distance of 2 compared photons falls within fuzzy distance range
					int distance = dist(photonY, photonX, centerY, centerX);
					if (distance >= radius * (1-fuzzy) && distance <= radius * (1+fuzzy))
						acc++;
				}
			}
			accumulator[centralPhoton] = acc;
		}

		printf("Accumulator completed: \n"); //outputs the 3 "most promising" annuli locations--for now, that's just the accumulator cells with the highest count
		int* halos = threelargest(accumulator, photonCount, radius, xs, ys);
		if (halos[0] > maxAcc) {
			maxAcc = halos[0];
			argmaxX = halos[1];
			argmaxY = halos[2];
			argmaxRad = radius;
		}
		free(halos);
		free(accumulator);
	}

	int *ret = malloc(3 * sizeof(int));
	ret[0] = argmaxX;
	ret[1] = argmaxY;
	ret[2] = argmaxRad;
	return ret;
}
