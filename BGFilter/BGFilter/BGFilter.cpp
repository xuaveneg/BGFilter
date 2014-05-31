#include <cmath>
#include <iostream>

#define GAUSS(kernelCarre, x, y) \
	(float)((1 / (2 * 3.14159265358979323846 * (kernelCarre))) * exp(-((x) * (x) + (y) * (y)) / (2 * (kernelCarre))))

#define RET_OK				0
#define RET_KO_INPUTN		001
#define RET_KO_INPUTTYPE	002
#define RET_KO_INPUTFILE	003

#define RET_KO_MATRIX_INCOMPSIZE	999
#define RET_KO_MATRIX_OOB			999


class Matrix {
private:
	float *_m;
public:
	const unsigned int w;
	const unsigned int h;
	float &operator()(unsigned int i, unsigned int j) {
		if (i > w || j > h)
			throw RET_KO_MATRIX_OOB;
		return _m[j*w + i];
	}
	const float &operator()(unsigned int i, unsigned int j) const {
		if (i > w || j > h)
			throw RET_KO_MATRIX_OOB;
		return _m[j*w + i];
	}
	const Matrix& operator=(const Matrix& assign) {
		if (w != assign.w || h != assign.h)
			throw RET_KO_MATRIX_INCOMPSIZE;
		memcpy(_m, assign._m, w*h*sizeof(float));
		return *this;
	}
	Matrix(const Matrix &copy) :
		w(copy.w),
		h(copy.h) {
		_m = new float[w*h];
		memcpy(_m, copy._m, w*h*sizeof(float));
	}
	Matrix() :
		w(0),
		h(0) {
		_m = new float[w*h];
		memset(_m, 0, w*h*sizeof(float));
	}
	Matrix(unsigned int width, unsigned int height) :
		w(width),
		h(height) {
		_m = new float[w*h];
		memset(_m, 0, w*h*sizeof(float));
	}
	Matrix(unsigned int size) :
		w(size),
		h(size) {
		_m = new float[w*h];
		memset(_m, 0, w*h*sizeof(float));
	}
	Matrix(float *matrix, unsigned int width, unsigned int height) :
		h(height),
		w(width) {
		_m = new float[w*h];
		memcpy(_m, matrix, w*h*sizeof(float));
	}
	Matrix(float *matrix, unsigned int size) :
		h(size),
		w(size) {
		_m = new float[w*h];
		memcpy(_m, matrix, w*h*sizeof(float));
	}
	~Matrix() {
		delete[] _m;
	}
};

class GBFilter {
public:
	GBFilter(float kernelSize, int tileWidth, int tileHeight) :
		_kernelSize(kernelSize),
		_tileWidth(tileWidth),
		_tileHeight(tileHeight) {
		int matrixSize = (int)(6.f * kernelSize);
		_convolutionMatrix = new Matrix(matrixSize);
		for (int j = 0; j < matrixSize; ++j) {
			for (int i = 0; i < matrixSize; ++i) {
				(*_convolutionMatrix)(i, j) = GAUSS(kernelSize * kernelSize, i - matrixSize / 2, j - matrixSize / 2);
			}
		}
	}
	const unsigned int convSize() const {
		return _convolutionMatrix->w;
	}
	const float operator()(int i, int j) const {
		return (*_convolutionMatrix)(i, j);
	}
	const unsigned int tileWidth() const {
		return _tileWidth;
	}
	const unsigned int tileHeight() const {
		return _tileHeight;
	}
private:
	const float _kernelSize;
	const unsigned int _tileWidth;
	const unsigned int _tileHeight;
	Matrix *_convolutionMatrix;

	//friend class Bitmap;
// We don't want GBFilter to be initialized from nothing nor copied
private:
	GBFilter() : _kernelSize(0), _tileWidth(0), _tileHeight(0){};
	GBFilter(const GBFilter &filter) : _kernelSize(filter._kernelSize), _tileWidth(filter._tileWidth), _tileHeight(filter._tileHeight) {};
	const GBFilter &operator=(const GBFilter &filter) {}
};

class Bitmap {
public:
	const unsigned int w;
	const unsigned int h;
	Bitmap(const char *inputFile) :
		w(0),
		h(0) {
		char c;
		FILE *f = fopen(inputFile, "r+");
		while (0 != fread(&c, 1, 1, f)) {
			std::cout << (int)c;
		}
	}
	void apply(const GBFilter &filter) {
		unsigned int convSize = filter.convSize();
		unsigned int iOut, jOut;
		unsigned int tileWidth = filter.tileWidth(), tileHeight = filter.tileHeight();
		for (unsigned int color = 0; color < 3; ++color) {
			Matrix outMatrix(w, h);
			Matrix inMatrix(w, h);// *(data[color]));
			for (unsigned int jIn = 0; jIn < h; ++jIn) {
				for (unsigned int iIn = 0; iIn < w; ++iIn) {
					for (unsigned int jConv = 0; jConv < convSize; ++jConv) {
						jOut = jIn - (jConv - convSize / 2);
						if (jOut / tileHeight != jIn / tileHeight)
							// The column is out of computing bound for this tile (note : works also for < 0 with overflow)
							break;
						for (unsigned int iConv = 0; iConv < convSize; ++iConv) {
							iOut = iIn - (iConv - convSize / 2);
							if (iOut / tileWidth != iIn / tileWidth)
								// The row is out of computing bound for this tile (note : works also for < 0 with overflow)
								break;
							outMatrix(0, 0) += inMatrix(iIn, jIn) * filter(iConv, jConv);
						}
					}
				}
			}
			//(*data[color]) = outMatrix;
		}
	}
	void write(const char* outputFile) const {

	}
private:
	Matrix *data[3];
// We don't want Bitmap to be initialized from nothing nor copied
private:
		Bitmap() : w(0), h(0) {}
		Bitmap(const Bitmap &b) : w(b.w), h(b.h) {}
		const Bitmap &operator=(const Bitmap &bitmap) {}

};

int main(int argc, char **argv)
{
	char *input_file;
	char *output_file;
	float kernel_size;
	int tile_width;
	int tile_height;
	if (argc == 6) {
		input_file = argv[1];
		output_file = argv[2];
		kernel_size = (float)atof(argv[3]);
		tile_width = atoi(argv[4]);
		tile_height = atoi(argv[5]);
	}
	else if (argc == 4) {
		input_file = argv[1];
		output_file = argv[2];
		kernel_size = (float)atof(argv[3]);
	}
	else {
		std::cout << "Usage:" << std::endl;
		std::cout << "gbfilter input_file output_file kernel_size tile_width tile_height" << std::endl;
		std::cout << std::endl;
		std::cout << "Example:" << std::endl;
		std::cout << "gbfilter my_image.bmp my_result.bmp 50.2 64 64" << std::endl;
		std::cout << std::endl;
		std::cout << "Description:" << std::endl;
		std::cout << "input_file" << std::endl;
		std::cout << "24 - bit BMP input image file." << std::endl;
		std::cout << "output_file" << std::endl;
		std::cout << "24 - bit BMP output image file." << std::endl;
		std::cout << "kernel_size" << std::endl;
		std::cout << "set in pixel or sub - pixel the size of the kernel." << std::endl;
		std::cout << "tile_width" << std::endl;
		std::cout << "set the width of the tile" << std::endl;
		std::cout << "tile_height" << std::endl;
		std::cout << "set the height of the tile" << std::endl;
		return RET_KO_INPUTN;
	}
	GBFilter filter(kernel_size, tile_width, tile_height);
	Bitmap bitmap(input_file);
	bitmap.apply(filter);
	bitmap.write(output_file);
	return RET_OK;
}