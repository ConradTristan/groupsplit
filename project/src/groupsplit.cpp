#include <iostream>
#include <fstream>
#include <ctype.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <array>
#include <new> 
#include <memory>
#include <algorithm>

using namespace std;

struct Pixel {
	unsigned char red, green, blue;
	float sum;
};

struct PixelMap {
	bool* pixels;
};

struct Group {
	float average;
	float min;
	int numPixels;
	bool* pixelMap;
	bool isPartOfChain;
};

struct Block {
	Pixel* pixels;
	int pixelsLength;
	Group* groups;
};

class RowNode {
    public:
        Group group;
    RowNode* next;
    int totalX;
};

class ColumnNode {
    public:
        RowNode* startNode;
    ColumnNode* next;
    int totalY;
};

class GroupNode {
    public:
        ColumnNode* startNode;
    GroupNode* next;
    int totalPixels;
    float average;
    float min;
};

// global variables to prevent pointers addressing garbage after functions

streampos size;
int colorDataIndex;

char* vals;
Pixel* pixels;
Block* blocks;
GroupNode* groupNodeList;

PixelMap* pixelMaps;

int width, height, maxval;

int blockWidth, blockHeight;
int finalColumnWidth, finalRowHeight;

float maxSumValue = 0, minSumValue = 765;

// arguments 
char* filename;
char* outputname;
// with default values
int blockSize = 16, maxGroups = 4, groupSumThreshold = 50, maxColors = 2;

void read_ppm(string filename) {
	// ios::ate positions the get pointer at EOF
	ifstream ppmFile(filename, ios::in | ios::binary | ios::ate);
	if (ppmFile.is_open()) {
		// pointer is at EOF, so size of whole file is equal to pointer index value
		size = ppmFile.tellg();	
		vals = new char[(int) size];
		// reads entire file into vals
		ppmFile.seekg(0, ios::beg);
		ppmFile.read(vals, size);
		ppmFile.close();
	} else {
		cout << "Failed to read file!" << endl;
		exit(1);
	}
}

// checks if a byte is a digit under ascii interpretation
bool is_ascii_number(char n) {
	return n >= '0' && n <= '9';
}

void parse_ppm() {
	// check for valid magic number P6
	if (vals[0] == 'P' && vals[1] == '6') {
		int index = 2;
		// traverse whitespace to next block
		while (isspace(vals[index])) {
			index++;
		}
		// width definition (as ASCII)
		width = 0;
		while (!isspace(vals[index])) {	
			width *= 10;
			if (is_ascii_number(vals[index])) {
				width += vals[index] - 48; // ascii to number conversion
				index++;
			} else {
				cout << "Invalid PPM format!" << endl;
				exit(1);
			}
		}
		// traverse whitespace to next block
		while (isspace(vals[index])) {
			index++;
		}
		// height definition (as ASCII)
		height = 0;
		while (!isspace(vals[index])) {
			height *= 10;
			if (is_ascii_number(vals[index])) {
				height += vals[index] - 48; // ascii to number conversion
				index++;
			} else {
				cout << "Invalid PPM format" << endl;
				exit(1);
			}
		}
		// traverse whitespace to next block
		while (isspace(vals[index])) {
			index++;
		}
		// height definition (as ASCII)
		maxval = 0;
		while (!isspace(vals[index])) {
			maxval *= 10;
			if (is_ascii_number(vals[index])) {
				maxval += vals[index] - 48; // ascii to number conversion
				index++;
			} else {
				cout << "Invalid PPM format" << endl;
				exit(1);
			}
		}
		// traverse whitespace to next block
		while (isspace(vals[index])) {
			index++;
		}
		colorDataIndex = index;
		// traverse color data and store it in pixels array
		pixels = new Pixel[width * height];
		int pixelIndex = 0;
		while (index < (int) size) {
			unsigned char red = vals[index];
			index++;
			unsigned char green = vals[index];
			index++;
			unsigned char blue = vals[index];
			float sum = red + green + blue;
			Pixel pix = {
				red,
				green,
				blue,
				sum,
			};
			if (sum < minSumValue) {
				minSumValue = sum;
			}
			if (sum > maxSumValue) {
				maxSumValue = sum;
			}
			pixels[pixelIndex] = pix;
			pixelIndex++;
			index++;
		}
	} else {
		cout << "Invalid PPM format!" << endl;
		exit(1);
	}
}

void blockify_pixels() {
	
	// calculate how many blocks we subdivide the image into
	blockWidth = width / blockSize;
	blockWidth += ((width % blockSize == 0) ? 0 : 1);
	blockHeight = height / blockSize;
	blockHeight += ((height % blockSize == 0) ? 0 : 1);
	
	// the final row/column may be smaller, so we store those values for index accessing
	finalColumnWidth = width % blockSize;
	if (finalColumnWidth == 0) {
		finalColumnWidth = blockSize;
	}
	
	finalRowHeight = height % blockSize;
	if (finalRowHeight == 0) {
		finalRowHeight = blockSize;
	}
	
	blocks = new Block[blockWidth * blockHeight];
	
	// initialise blocks and add the appropriate pixels to its array
#pragma omp parallel
{
#pragma omp for
	for (int i = 0; i < blockWidth * blockHeight; i++) {
		
		int totalX = i % blockWidth;
		int totalY = i / blockWidth;
		
		int maxBlockX = blockSize;
		int maxBlockY = blockSize;
		
		if (totalX == blockWidth - 1) {
			maxBlockX = finalColumnWidth;
		} 
		if (totalY == blockHeight - 1) {	
			maxBlockY = finalRowHeight;
		}
		
		Block block = {
			new Pixel[maxBlockX*maxBlockY],
			maxBlockX*maxBlockY,
			new Group[maxGroups],
		};
		for (int blockY = 0; blockY < maxBlockY; blockY++) {
			for (int blockX = 0; blockX < maxBlockX; blockX++) {
				block.pixels[blockY*maxBlockX + blockX] = pixels[totalY * width * blockSize
				                       + totalX * blockSize
				                       + blockY * width
				                       + blockX];
			}
		}
		blocks[totalY*blockWidth + totalX] = block;	
	}
}
}

// adds a pixel to a group, calculating new average/minimum and updating the pixelmap
void update_group_with_pixel(Group &group, Pixel &pixel, int pixelIndex) {
	group.average *= group.numPixels;
	group.min = min(group.min, pixel.sum);
	group.numPixels++;
	group.average += pixel.sum;
	group.average /= group.numPixels;
	group.pixelMap[pixelIndex] = true;
}

void assign_groups_in_blocks() {
	pixelMaps = new PixelMap[blockWidth * blockHeight * maxGroups];
	// iterate over each block
#pragma omp parallel
{
#pragma omp for
	for (int i = 0; i < blockWidth * blockHeight; i++) {	
		// initalise each group with negative average and empty pixelMap
		// initalise it as part of a chain so we skip it later if it never receives a pixel
		for (int j = 0; j < maxGroups; j++) {
			PixelMap pixelMap = {
				new bool[blocks[i].pixelsLength] {},
			};
			pixelMaps[i*j+j] = pixelMap;
			Group group = {
				-1.0,
				755.0,
				0,	
				pixelMaps[i*j+j].pixels,
				true,
			};
			blocks[i].groups[j] = group;
		}
		// iterate over all pixels in a block to assign them to groups
		for (int j = 0; j < blocks[i].pixelsLength; j++) {
			bool groupAssigned = 0;
			int closestGroupIndex = 0;
			float diffToClosestGroup = 100000.0;
			// iterate over all groups
			for (int k = 0; k < maxGroups; k++) {
				float diffFromAverage = abs(blocks[i].groups[k].average - blocks[i].pixels[j].sum);
				// if we come across a group without other pixels, initalise it with the current pixel
				if (blocks[i].groups[k].average < 0) {
					blocks[i].groups[k].average = blocks[i].pixels[j].sum;
					blocks[i].groups[k].min = blocks[i].pixels[j].sum;
					blocks[i].groups[k].numPixels++;
					blocks[i].groups[k].pixelMap[j] = true;
					blocks[i].groups[k].isPartOfChain = false;
					groupAssigned = 1;
					break;
				// if a pixel is within the threshold of a group, add it there
				} else if (diffFromAverage < groupSumThreshold) {
					update_group_with_pixel(blocks[i].groups[k], blocks[i].pixels[j], j);	
					groupAssigned = 1;
					break;
				// if the group was not within the threshold, keep track of whether it was the closest group so far
				} else if (diffFromAverage < diffToClosestGroup) {;
					diffToClosestGroup = diffFromAverage;
					closestGroupIndex = k;
				}
			}
			// if we have filled all groups and none is within the threshold, add the pixel to the group with closest value
			if (!groupAssigned) {
				update_group_with_pixel(blocks[i].groups[closestGroupIndex], blocks[i].pixels[j], j);
				groupAssigned = 1;
			}
		}
	}
}	
}

void merge_groups() {
	groupNodeList = NULL;
	
	GroupNode* currentGroupNode;
	currentGroupNode = new GroupNode();
	// iterate over all blocks
	for (int startBlockIndex = 0; startBlockIndex < blockWidth * blockHeight; startBlockIndex++) {	
		// iterate over all groups in each block
		for (int j = 0; j < maxGroups; j++) {
			// if current group currently isn't in a chain, start a new chain
			if (!blocks[startBlockIndex].groups[j].isPartOfChain) {
				
				RowNode* currentNode;				
				currentNode = new RowNode();
				currentNode -> group = blocks[startBlockIndex].groups[j];
				currentNode -> totalX = startBlockIndex % blockWidth;
				blocks[startBlockIndex].groups[j].isPartOfChain = true;
				currentNode -> next = NULL;
				
				ColumnNode* currentRowStartNode;
				currentRowStartNode = new ColumnNode();
				currentRowStartNode -> startNode = currentNode;
				currentRowStartNode -> totalY = startBlockIndex / blockWidth;
				currentRowStartNode -> next = NULL;
				
				GroupNode* newGroupNode;
				newGroupNode = new GroupNode();
				newGroupNode -> startNode = currentRowStartNode;
				newGroupNode -> next = NULL;
				newGroupNode -> totalPixels = currentNode -> group.numPixels;
				newGroupNode -> average = currentNode -> group.average;
				newGroupNode -> min = currentNode -> group.min;
				
				// initalise final group list if not yet initialised
				if (groupNodeList == NULL) {
					currentGroupNode = newGroupNode;
					groupNodeList = currentGroupNode;
				// otherwise append it to the list
				} else {
					currentGroupNode -> next = newGroupNode;
					currentGroupNode = newGroupNode;
				}				
				
				// iterate over all rows starting at blockIndex		
				for (int totalY = startBlockIndex / blockWidth; totalY < blockHeight; totalY++) {
					currentNode = currentRowStartNode -> startNode;
					// if we aren't in the starting row, see if the starting block in the next row can link
					if (!(totalY == startBlockIndex / blockWidth)) {
						// iterate over all groups in the block in the next row
						for (int currentGroupIndex = 0; currentGroupIndex < maxGroups; currentGroupIndex++) {
							// if a group is close enough, link it to the group and continue from that block
							if (abs(currentNode -> group.average - blocks[totalY*blockWidth + startBlockIndex % blockWidth].groups[currentGroupIndex].average) < groupSumThreshold) {
								if (!blocks[totalY*blockWidth + startBlockIndex % blockWidth].groups[currentGroupIndex].isPartOfChain) {
									
									RowNode* newNode;
									newNode = new RowNode();
									newNode -> group = blocks[totalY*blockWidth + startBlockIndex % blockWidth].groups[currentGroupIndex];
									newNode -> totalX = startBlockIndex % blockWidth;
									blocks[totalY*blockWidth + startBlockIndex % blockWidth].groups[currentGroupIndex].isPartOfChain = true;
									newNode -> next = NULL;;
									currentNode = newNode;
																	
									currentGroupNode -> average *= currentGroupNode -> totalPixels;
									currentGroupNode -> totalPixels += currentNode -> group.numPixels;
									currentGroupNode -> average += currentNode -> group.average * currentNode -> group.numPixels;
									currentGroupNode -> average /= currentGroupNode -> totalPixels;
									currentGroupNode -> min = min(currentGroupNode -> min, newNode -> group.min);
									
									ColumnNode* newColumnNode;
									newColumnNode = new ColumnNode();
									newColumnNode -> startNode = currentNode;
									newColumnNode -> totalY = totalY;
									newColumnNode -> next = NULL;
									
									currentRowStartNode -> next = newColumnNode;
									currentRowStartNode = newColumnNode;							
									
									goto iterateX;
								}
							}
						}
						// if it can't link, the linking of this group is finished
						goto complete;
					}
					iterateX:;
					// iterate over all following blocks in current row
					for (int totalX = startBlockIndex % blockWidth + 1; totalX < blockWidth; totalX++) {
						// iterate over all groups in next block
						for (int currentGroupIndex = 0; currentGroupIndex < maxGroups; currentGroupIndex++) {
							// if a group is close enough, link it to the group and continue from that block
							if (abs(currentNode -> group.average - blocks[totalY*blockWidth + totalX].groups[currentGroupIndex].average) < groupSumThreshold) {
								if (!blocks[totalY*blockWidth + totalX].groups[currentGroupIndex].isPartOfChain) {
																		
									RowNode* newNode;
									newNode = new RowNode();
									newNode -> group = blocks[totalY*blockWidth + totalX].groups[currentGroupIndex];
									newNode -> totalX = totalX;
									blocks[totalY*blockWidth + totalX].groups[currentGroupIndex].isPartOfChain = true;
									newNode -> next = NULL;
									currentNode -> next = newNode;
									currentNode = newNode;
									
									currentGroupNode -> average *= currentGroupNode -> totalPixels;
									currentGroupNode -> totalPixels += currentNode -> group.numPixels;
									currentGroupNode -> average += currentNode -> group.average * currentNode -> group.numPixels;
									currentGroupNode -> average /= currentGroupNode -> totalPixels;
									currentGroupNode -> min = min(currentGroupNode -> min, newNode -> group.min);
																																		
									goto nextXBlock;
								}
							}
						}
						// if we don't find a group close enough, continue with the next row
						goto nextYBlock;
						nextXBlock:;
					}
					
					nextYBlock:;
				}			
				complete:;
			}
		}
	}					
}
		
void recolor() {
	
	// divide by 3 to bring them back into [0, 255]
	minSumValue /= 3;
	maxSumValue /= 3;
	GroupNode* currentGroupNode;
	ColumnNode* currentColumnNode;
	RowNode* currentRowNode;
	
	// iterate through all the linked groups to recolor their pixels
	currentGroupNode = groupNodeList;
	while (currentGroupNode != NULL) {
		currentColumnNode = currentGroupNode -> startNode;
		while (currentColumnNode != NULL) {
			currentRowNode = currentColumnNode -> startNode;
			while (currentRowNode != NULL) {
				int maxBlockX = blockSize;
				int maxBlockY = blockSize;
				
				if (currentRowNode -> totalX == blockWidth - 1) {
					maxBlockX = finalColumnWidth;
				} 
				if (currentColumnNode -> totalY == blockHeight - 1) {	
					maxBlockY = finalRowHeight;
				}
								
				float currentValue = currentGroupNode -> min;
				currentValue /= 3;
				
				// rescale values so that they more accurately go from white to black
				currentValue = (255 / (maxSumValue - minSumValue)) * (currentValue - minSumValue);
				int sectorLength = (maxval / (maxColors - 1));
				
				// Index of the section of color division in which the current pixel falls
				int recolorSector = currentValue / sectorLength;
				// Assign the index depending on which is closer in the current section
				int recolorIndex = (currentValue - recolorSector * sectorLength > sectorLength / 2) ? recolorSector + 1 : recolorSector;
				currentValue = recolorIndex * sectorLength;
				
				for (int blockY = 0; blockY < maxBlockY; blockY++) {
					for (int blockX = 0; blockX < maxBlockX; blockX++) {
						int pixelIndex = currentColumnNode -> totalY * width * blockSize
											   + currentRowNode -> totalX * blockSize
											   + blockY * width
											   + blockX;									        
				        
				        // Assign the color of all pixels in the current group to the new value
						if (currentRowNode -> group.pixelMap[blockY*maxBlockX + blockX]) {																				
							pixels[pixelIndex].red = currentValue;
							pixels[pixelIndex].green = currentValue;
							pixels[pixelIndex].blue = currentValue;
						}
					}
				}
				
				currentRowNode = currentRowNode -> next;
			}
			currentColumnNode = currentColumnNode -> next;
		}
		currentGroupNode = currentGroupNode -> next;
	}					
}

void write_ppm_file(string filename) {
	// write the new pixel data into vals
	int pixelIndex = 0;
	for (int i = colorDataIndex; i < size; i += 3) {
		vals[i] = pixels[pixelIndex].red;
		vals[i+1] = pixels[pixelIndex].green;
		vals[i+2] = pixels[pixelIndex].blue;
		pixelIndex++;
	}
	
	// write the full vals array to the outpute file
	ofstream ppmFileOut(filename, ios::out | ios::binary);
	if (ppmFileOut.is_open()) {
		ppmFileOut.write(vals, size);
		ppmFileOut.close();
	} else {
		cout << "Error writing file" << endl;
		exit(1);
	}
}

char* getCmdOption(char ** begin, char ** end, const std::string & option) {
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

int main(int argc, char* argv[]) {
	
	char* fileArg = getCmdOption(argv, argv + argc, "-f");
	char* outputArg = getCmdOption(argv, argv + argc, "-o");
	
	char* sizeArg = getCmdOption(argv, argv + argc, "-s");
	char* groupsArg = getCmdOption(argv, argv + argc, "-g");
	char* thresholdArg = getCmdOption(argv, argv + argc, "-t");
	char* colorsArg = getCmdOption(argv, argv + argc, "-c");
	
	if (fileArg) {
		filename = fileArg;
	} else {
		cout << "You must specify an input file!" << endl;
		exit(1);
	}
	
	if (outputArg) {
		outputname = outputArg;
	} else {
		cout << "You must specify an output file!" << endl;
		exit(1);
	}
	
	if (sizeArg) {
		blockSize = stoi(sizeArg);
	}
	
	if (groupsArg) {
		maxGroups = stoi(groupsArg);
	}
	
	if (thresholdArg) {
		groupSumThreshold = stoi(thresholdArg);
	}
	
	if (colorsArg) {
		maxColors = stoi(colorsArg);
	}
	
	double start = omp_get_wtime();
		
	read_ppm(filename);
	double read_time = omp_get_wtime() - start;
		
	parse_ppm();
	double parse_time = omp_get_wtime() - start;
		
	blockify_pixels();
	double blockify_time = omp_get_wtime() - start;
		
	assign_groups_in_blocks();
	double assign_time = omp_get_wtime() - start;

	merge_groups();
	double merge_time = omp_get_wtime() - start;

	recolor();
	double recolor_time = omp_get_wtime() - start;

	write_ppm_file(outputname);
	double write_time = omp_get_wtime() - start;
		
	cout << "Read: " << read_time << endl;
	cout << "Parse: " << parse_time << endl;
	cout << "Blockify: " << blockify_time << endl;
	cout << "Assign: " << assign_time << endl;
	cout << "Merge: " << merge_time << endl;
	cout << "Recolor: " << recolor_time << endl;
	cout << "Write: " << write_time << endl;
	
	return 0;
}
