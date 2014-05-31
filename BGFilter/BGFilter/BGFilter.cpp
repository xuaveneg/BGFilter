#include <cmath>
#include <cstdint>
#include <iostream>

#define GAUSS(kernelCarre, x, y) \
	(float)((1 / (2 * 3.14159265358979323846 * (kernelCarre))) * exp(-((x) * (x) + (y) * (y)) / (2 * (kernelCarre))))

#define RET_OK				0
#define RET_KO_INPUTN		001
#define RET_KO_INPUTTYPE	002
#define RET_KO_INPUTFILE	003

#define RET_KO_MATRIX_OOB			101

#define RET_KO_FILEREAD				201
#define RET_KO_FILEWRITE			202
#define RET_KO_FILEHEAD				203
#define RET_KO_FILEBODY				204
#define RET_KO_FILEUNSUPPORTED		205

// Matrix Utils
class Matrix {
private:
	float *_m; // Data
	unsigned int _w; // Width
	unsigned int _h; // Height
public:
	// Width and Height accessors
	inline const unsigned int w() const {
		return _w;
	}
	inline const unsigned int h() const {
		return _h;
	}
	// Data modifier
	float &operator()(unsigned int i, unsigned int j) {
		if (i > _w || j > _h)
			throw RET_KO_MATRIX_OOB;
		return _m[j*_w + i];
	}
	// Data accessor
	const float &operator()(unsigned int i, unsigned int j) const {
		if (i > _w || j > _h)
			throw RET_KO_MATRIX_OOB;
		return _m[j*_w + i];
	}
	// Assignment
	const Matrix& operator=(const Matrix& assign) {
		if (_w != assign._w || _h != assign._h) {
			_w = assign._w;
			_h = assign._h;
			delete[] _m;
			_m = new float[_w*_h];
		}
		memcpy(_m, assign._m, _w*_h*sizeof(float));
		return *this;
	}
	// Copy Constructor
	Matrix(const Matrix &copy) :
		_w(copy._w),
		_h(copy._h) {
		_m = new float[_w*_h];
		memcpy(_m, copy._m, _w*_h*sizeof(float));
	}
	// Default constructor : empty matrix
	Matrix() :
		_w(0),
		_h(0) {
		_m = new float[_w*_h];
		memset(_m, 0, _w*_h*sizeof(float));
	}
	// Constructor given its size : 0 matrix
	Matrix(unsigned int width, unsigned int height) :
		_w(width),
		_h(height) {
		_m = new float[_w*_h];
		memset(_m, 0, _w*_h*sizeof(float));
	}
	// Constructor given its size for square : 0 matrix
	Matrix(unsigned int size) :
		_w(size),
		_h(size) {
		_m = new float[_w*_h];
		memset(_m, 0, _w*_h*sizeof(float));
	}
	// Constructor given its size and data (no data check)
	Matrix(float *matrix, unsigned int width, unsigned int height) :
		_h(height),
		_w(width) {
		_m = new float[_w*_h];
		memcpy(_m, matrix, _w*_h*sizeof(float));
	}
	// Constructor given its size and data for square (no data check)
	Matrix(float *matrix, unsigned int size) :
		_h(size),
		_w(size) {
		_m = new float[_w*_h];
		memcpy(_m, matrix, _w*_h*sizeof(float));
	}
	// Destructor
	~Matrix() {
		delete[] _m;
	}
};

// Class describing Gaussian Blur Filter
class GBFilter {
public:
	// Construct a Gaussian Blur Filter : needs informations (we don't want to redefine it after initialization)
	GBFilter(float kernelSize, int tileWidth, int tileHeight) :
		_kernelSize(kernelSize),
		_tileWidth(tileWidth),
		_tileHeight(tileHeight) {
		std::cout << "Computing Convolution Matrix ... ";
		// Construct convolution matrix
		int matrixSize = (int)(6.f * kernelSize);
		_convolutionMatrix = new Matrix(matrixSize);
		for (int j = 0; j < matrixSize; ++j) {
			for (int i = 0; i < matrixSize; ++i) {
				(*_convolutionMatrix)(i, j) = GAUSS(kernelSize * kernelSize, i - matrixSize / 2, j - matrixSize / 2);
			}
		}
		std::cout << "Done." << std::endl;
		std::cout << "Convolution Matrix Size: " << matrixSize << std::endl;
	}
	// Destruct a Gaussian Blur Filter, ie its Convolution Matrix
	~GBFilter() {
		delete _convolutionMatrix;
	}
	// Convolution Informations accessors
	const unsigned int convSize() const {
		return _convolutionMatrix->w();
	}
	const float operator()(int i, int j) const {
		return (*_convolutionMatrix)(i, j);
	}
	// Filter informations accessors
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

// class describing Bitmap File
class Bitmap {
private:
	// File reading operations change with endianness
	bool bigEndian;
	// File reading utils
	// Reads an unsigned 32 bit value
	uint32_t fileRead32b(FILE *f) {
		uint32_t ret(0);
		char c;
		for (unsigned int i = 0; i < 32; i += 8) {
			if (1 != fread(&c, 1, 1, f)) throw RET_KO_FILEREAD;
			// As it is unsigned and as we know the file endianness, we don't care about endianness, just add the value
			ret += c << i;
		}
		return ret;
	}
	// Reads an unsigned 16 bit value
	uint16_t fileRead16b(FILE *f) {
		uint16_t ret(0);
		char c;
		for (unsigned int i = 0; i < 16; i += 8) {
			if (1 != fread(&c, 1, 1, f)) throw RET_KO_FILEREAD;
			// As it is unsigned and as we know the file endianness, we don't care about endianness, just add the value
			ret += c << i;
		}
		return ret;
	}
	// Reads a signed 32 bit value
	int32_t fileRead32bS(FILE *f) {
		int32_t ret(0);
		char *ret_c = (char*)(&ret);
		char c;
		for (unsigned int i = 0; i < 4; ++i) {
			if (1 != fread(&c, 1, 1, f)) throw RET_KO_FILEREAD;
			// Here we need endianness. Bitmap file uses little endian
			 ret_c[(bigEndian ? 3-i : i)] = c;
		}
		return ret;
	}
	// Reads a signed 16 bit value
	int16_t fileRead16bS(FILE *f) {
		int16_t ret(0);
		char *ret_c = (char*)(&ret);
		char c;
		for (unsigned int i = 0; i < 2; ++i) {
			if (1 != fread(&c, 1, 1, f)) throw RET_KO_FILEREAD;
			// Here we need endianness. Bitmap file uses little endian
			ret_c[(bigEndian ? 1 - i : i)] = c;
		}
		return ret;
	}

	// Reads a pixel data
	void fileReadPixel(float &red, float &green, float &blue, FILE *f) {
		if (1 != fread(&blue, 1, 1, f)) throw RET_KO_FILEREAD;
		if (1 != fread(&green, 1, 1, f)) throw RET_KO_FILEREAD;
		if (1 != fread(&red, 1, 1, f)) throw RET_KO_FILEREAD;
	}
	// File writing utils
	// Writes an unsigned 32 bit value
	void fileWrite32b(uint32_t value, FILE *f) const {
		char c;
		for (unsigned int i = 0; i < 32; i += 8) {
			// As it is unsigned and as we know the file endianness, we don't care about endianness, just get the right byte
			c = value >> i;
			if (1 != fwrite(&c, 1, 1, f)) throw RET_KO_FILEWRITE;
		}
	}
	// Writes an unsigned 16 bit value
	void fileWrite16b(uint16_t value, FILE *f) const {
		char c;
		for (unsigned int i = 0; i < 16; i += 8) {
			// As it is unsigned and as we know the file endianness, we don't care about endianness, just get the right byte
			c = value >> i;
			if (1 != fwrite(&c, 1, 1, f)) throw RET_KO_FILEWRITE;
		}
	}
	// Writes a signed 32 bit value
	void  fileWrite32bS(int32_t value, FILE *f) const {
		char *val_c = (char*)(&value);
		char c;
		for (unsigned int i = 0; i < 4; ++i) {
			// Here we need endianness. Bitmap file uses little endian
			c = val_c[(bigEndian ? 3 - i : i)];
			if (1 != fwrite(&c, 1, 1, f)) throw RET_KO_FILEWRITE;
		}
	}
	// Writes a signed 16 bit value
	void fileWrite16bS(int16_t value, FILE *f) const {
		char *val_c = (char*)(&value);
		char c;
		for (unsigned int i = 0; i < 2; ++i) {
			// Here we need endianness. Bitmap file uses little endian
			c = val_c[(bigEndian ? 1 - i : i)];
			if (1 != fwrite(&c, 1, 1, f)) throw RET_KO_FILEWRITE;
		}
	}
	// Writes a pixel data
	void fileWritePixel(float red, float green, float blue, FILE *f) const {
		if (1 != fwrite(&blue, 1, 1, f)) throw RET_KO_FILEWRITE;
		if (1 != fwrite(&green, 1, 1, f)) throw RET_KO_FILEWRITE;
		if (1 != fwrite(&red, 1, 1, f)) throw RET_KO_FILEWRITE;
	}
	// Bitmap File structs
	struct {
		uint16_t magicNumber	= 0x4D42; // Bitmap File is little endian => "BM" (0x424D) is read as 0X4D42
		uint32_t fileSize;
		uint16_t reserved1;
		uint16_t reserved2;
		uint32_t dataOffset;
	} fileHeader;
	struct {
		uint32_t headerSize;
		int32_t width;
		int32_t height;
		uint16_t colorPlanes	= 1;
		uint16_t bitsPerPixel	= 24;
		uint32_t compression	= 0;
		uint32_t rawSize;
		int32_t horizontalRes;
		int32_t verticalRes;
		uint32_t nColor			= 0;
		uint32_t nColorUsed		= 0;
	} informationHeader;
	void *additionalHeader;
	uint32_t _w;
	uint32_t _h;
public:
	Bitmap(const char *inputFile) {
		std::cout << "Reading File ... ";
		// Read 24-bit Bitmap File
		FILE *f = fopen(inputFile, "rb");
		// Thanks, Wikipedia !
		// Bitmap File Header
		//fileReadHeadPattern("BM", 2, f);
		fileHeader.magicNumber  = fileRead16b(f);		// Identify BMP file	(offset 0)	(size 2)
		fileHeader.fileSize		= fileRead32b(f);		// BMP size				(offset 2)	(size 4)
		fileHeader.reserved1	= fileRead16b(f);		// Reserved				(offset 6)	(size 2)
		fileHeader.reserved2	= fileRead16b(f);		// Reserved				(offset 8)	(size 2)
		fileHeader.dataOffset	= fileRead32b(f);		// Pixel Array Offset	(offset 10)	(size 4)
		// Bitmap Information Header (supposed at least BITMAPINFOHEADER with Windows NT & 3.1x or later)
		informationHeader.headerSize	= fileRead32b(f);	// Size of this header (depend on version)	(offset 14)	(size 4)
		informationHeader.width			= fileRead32bS(f);	// Bitmap Width in pxls (signed int)		(offset 18) (size 4)
		informationHeader.height		= fileRead32bS(f);	// Bitmap Height in pxls (signed int)		(offset 22) (size 4)
		informationHeader.colorPlanes	= fileRead16b(f);	// Number of Color Planes (always 1)		(offset 26) (size 2)
		if (informationHeader.colorPlanes	!= 1)	throw RET_KO_FILEHEAD;
		informationHeader.bitsPerPixel	= fileRead16b(f);	// Number of bits per pxl (only 24)			(offset 28) (size 2)
		if (informationHeader.bitsPerPixel	!= 24)	throw RET_KO_FILEUNSUPPORTED;
		informationHeader.compression	= fileRead32b(f);	// Compression Method (only uncompressed)	(offset 30)	(size 4)
		if (informationHeader.compression	!= 0)	throw RET_KO_FILEUNSUPPORTED;
		informationHeader.rawSize		= fileRead32b(f);	// Raw bitmap data size						(offset 34)	(size 4)
		informationHeader.horizontalRes	= fileRead32bS(f);	// Horizontal Resolution					(offset 38)	(size 4)
		informationHeader.verticalRes	= fileRead32bS(f);	// Vertical Resolution						(offset 42)	(size 4)
		informationHeader.nColor		= fileRead32b(f);	// Color Palette Size (unsupported)			(offset 46)	(size 4)
		if (informationHeader.nColor != 0) throw RET_KO_FILEUNSUPPORTED;
		informationHeader.nColorUsed = fileRead32b(f);		// Important Colors Nukber (unsupported)	(offset 50) (size 4)
		if (informationHeader.nColorUsed != 0) throw RET_KO_FILEUNSUPPORTED;
		// Additional Headers (structure not known, kept as in) (offset 54)
		additionalHeader = malloc(fileHeader.dataOffset - 54);
		fread(additionalHeader, 1, fileHeader.dataOffset - 54, f);
		// Skip to Data
		_w = abs(informationHeader.width);
		_h = abs(informationHeader.height);
		long int tell;
		long int val = fileHeader.dataOffset;
		for (unsigned int i = 0; i < 3; ++i)
			data[i] = new Matrix(_w, _h);
		for (unsigned int j = 0; j < _h; ++j) {
			for (unsigned int i = 0; i < _w; ++i) {
				fileReadPixel((*data[0])(i, j), (*data[1])(i, j), (*data[2])(i, j), f);
				tell = ftell(f);
				val += 3;
			}
			tell = ftell(f);
			val += (4 - ((val - fileHeader.dataOffset) % 4)) % 4;
			fseek(f, (4 - ((tell - fileHeader.dataOffset) % 4)) % 4, SEEK_CUR);
		}
		fclose(f);
		std::cout << "Done." << std::endl;
		std::cout << "Size: " << _w << 'x' << _h << std::endl;
	}
	void apply(const GBFilter &filter) {
		std::cout << "Applying Filter ... " << std::endl;
		unsigned int convSize = filter.convSize();
		unsigned int iOut, jOut;
		unsigned int tileWidth = filter.tileWidth(), tileHeight = filter.tileHeight();
		for (unsigned int color = 0; color < 3; ++color) {
			std::cout << "Color " << color << " ... ";
			Matrix *outMatrix = new Matrix(_w, _h);
			Matrix *inMatrix = data[color];
			for (unsigned int jIn = 0; jIn < _h; ++jIn) {
				for (unsigned int iIn = 0; iIn < _w; ++iIn) {
					for (unsigned int jConv = 0; jConv < convSize; ++jConv) {
						jOut = jIn - (jConv - convSize / 2);
						if (jOut / tileHeight != jIn / tileHeight)
							// The column is out of computing bound for this tile (note : works also for < 0 with overflow)
							continue;
						for (unsigned int iConv = 0; iConv < convSize; ++iConv) {
							iOut = iIn - (iConv - convSize / 2);
							if (iOut / tileWidth != iIn / tileWidth)
								// The row is out of computing bound for this tile (note : works also for < 0 with overflow)
								continue;
							(*outMatrix)(0, 0) += (*inMatrix)(iIn, jIn) * filter(iConv, jConv);
						}
					}
				}
			}
			data[color] = outMatrix;
			delete inMatrix;
			std::cout << "Done." << std::endl;
		}
		std::cout << "Filter Applied." << std::endl;
	}
	void write(const char* outputFile) const {
		std::cout << "Writing Output File ... ";
		FILE *f = fopen(outputFile, "wb");
		// Bitmap File Header
		fileWrite16b(fileHeader.magicNumber, f);// Identify BMP file	(offset 0)	(size 2)
		fileWrite32b(fileHeader.fileSize, f);	// BMP size				(offset 2)	(size 4)
		fileWrite16b(fileHeader.reserved1, f);	// Reserved				(offset 6)	(size 2)
		fileWrite16b(fileHeader.reserved2, f);	// Reserved				(offset 8)	(size 2)
		fileWrite32b(fileHeader.dataOffset, f);	// Pixel Array Offset	(offset 10)	(size 4)
		// Bitmap Information Header (supposed at least BITMAPINFOHEADER with Windows NT & 3.1x or later)
		fileWrite32b(informationHeader.headerSize, f);		// Size of this header (depend on version)	(offset 14)	(size 4)
		fileWrite32bS(informationHeader.width, f);			// Bitmap Width in pxls (signed int)		(offset 18) (size 4)
		fileWrite32bS(informationHeader.height, f);			// Bitmap Height in pxls (signed int)		(offset 22) (size 4)
		fileWrite16b(informationHeader.colorPlanes, f);		// Number of Color Planes (always 1)		(offset 26) (size 2)
		fileWrite16b(informationHeader.bitsPerPixel, f);	// Number of bits per pxl (only 24)			(offset 28) (size 2)
		fileWrite32b(informationHeader.compression, f);		// Compression Method (only uncompressed)	(offset 30)	(size 4)
		fileWrite32b(informationHeader.rawSize, f);			// Raw bitmap data size						(offset 34)	(size 4)
		fileWrite32bS(informationHeader.horizontalRes, f);	// Horizontal Resolution					(offset 38)	(size 4)
		fileWrite32bS(informationHeader.verticalRes, f);	// Vertical Resolution						(offset 42)	(size 4)
		fileWrite32b(informationHeader.nColor, f);			// Color Palette Size (unsupported)			(offset 46)	(size 4)
		fileWrite32b(informationHeader.nColorUsed, f);		// Important Colors Nukber (unsupported)	(offset 50) (size 4)
		// Additional Headers (structure not known, kept as in) (offset 54)
		fwrite(additionalHeader, fileHeader.dataOffset - 54, 1, f);
		long int tell;
		const char *padding = "\0\0\0\0";
		for (unsigned int j = 0; j < _h; ++j) {
			for (unsigned int i = 0; i < _w; ++i) {
				fileWritePixel((*data[0])(i, j), (*data[1])(i, j), (*data[2])(i, j), f);
			}
			tell = ftell(f);
			fwrite(padding, 1, (4 - ((tell - fileHeader.dataOffset) % 4)) % 4, f);
		}
		if (fileHeader.fileSize != ftell(f)) throw RET_KO_FILEWRITE;
		fclose(f);
		std::cout << "Done." << std::endl;
	}
private:
	Matrix *data[3];
// We don't want Bitmap to be initialized from nothing nor copied
private:
	Bitmap() {}
	Bitmap(const Bitmap &b) {}
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