#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>

#define pi 3.1415926
#define JPEG_BLOCK_SIZE 8 // no change for current
using namespace std;

void mulmat(double A[], double mat[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE], double B[]) // transform array by multiple to matrix
{
    double temp = 0;
    for (int i=0; i<JPEG_BLOCK_SIZE; i++)
    {
        temp = 0;
        for (int j=0; j<JPEG_BLOCK_SIZE; j++)
            temp += A[j] * mat[i][j];
        B[i] = temp;
    }
}

void dct(char origin[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE], double dct_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE]) // DCT transform
{
    double temp[JPEG_BLOCK_SIZE*JPEG_BLOCK_SIZE];
    double dct_mat[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE];
    for (int i=0; i<JPEG_BLOCK_SIZE; i++) // DC part
        dct_mat[0][i] = 1/sqrt(2.0)/2.0;
    for (int i=1; i<JPEG_BLOCK_SIZE; i++) // AC part
        for (int j=1; j<=JPEG_BLOCK_SIZE; j++)
            dct_mat[i][j-1] = cos(i*pi*(2*j-1)/(2*JPEG_BLOCK_SIZE))/2.0;

    // TD
    double td[JPEG_BLOCK_SIZE];
    double fd[JPEG_BLOCK_SIZE];
    // First
    for (int i=0; i<JPEG_BLOCK_SIZE; i++)
    {
        for (int j=0; j<JPEG_BLOCK_SIZE; j++)
            td[j] = (double)origin[i][j];
        mulmat(td, dct_mat, fd);
        for (int j=0; j<JPEG_BLOCK_SIZE; j++)
            temp[i*JPEG_BLOCK_SIZE+j] = fd[j];
    }

    // Second
    for (int j=0; j<JPEG_BLOCK_SIZE; j++)
    {
        for (int i=0; i<JPEG_BLOCK_SIZE; i++)
            td[i] = temp[i*JPEG_BLOCK_SIZE+j];
        mulmat(td, dct_mat, fd);
        for (int i=0; i<JPEG_BLOCK_SIZE; i++)
            dct_r[i][j] = fd[i];
    }
}

void quant(double dct_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE], int quant_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE]) // quant the DC transform result
{
    int Q_mat[8][8] = {{16,11,10,16,24,40,51,61},
                    {12,12,14,19,26,58,60,55},
                    {14,13,16,24,40,57,69,56},
                    {14,17,22,29,51,87,80,62},
                    {18,22,37,56,68,109,103,77},
                    {24,35,55,64,81,104,113,92},
                    {49,64,78,87,103,121,120,101},
                    {72,92,95,98,112,100,103,99}};
    double temp;
    for (int i=0; i<JPEG_BLOCK_SIZE; i++)
        for (int j=0; j<JPEG_BLOCK_SIZE; j++)
            quant_r[i][j] = dct_r[i][j]/Q_mat[i][j] + ((dct_r[i][j] > 0)?0.5:(dct_r[i][j] < 0)?(-0.5):0.0);
}

void dc_entro(char code[], int PRE_DC, int DC) // entropy code the DC component
{
    char dc_len_table[12][10] = {"00","010", "011","100","101","110","1110","11110","111110","1111110","11111110","111111110"};
    int diff = DC - PRE_DC; // diff the DC component with previous block
    int len;
    if (diff == 0) // 0 can't do log2 calculation
        len = 0;
    else
        len = (int)floor(log2(abs(diff))) + 1; // length of binary of absolute value
    strcat(code, dc_len_table[len]);
    int unit = (diff>0)?diff:(diff+(2<<(len-1))-1); // mapping (-[2^(len)], -[2^(len-1)]) to (0, 2^(len-1)-1)
    char temp[16];
    memset(temp, 0, 16*sizeof(char));
    //for (int tmp = unit, i = 0; tmp >= 0; tmp /= 2, i++)
    //    temp[i] = (tmp%2)?'1':'0';
    //for (int tmp = unit, i = len-1; i>=0; tmp /= 2, i--)
    //    temp[i] = (tmp%2)?'1':'0';
    //temp[len] = '\0';
    for (int tmp = unit, i = len-1; i>=0; tmp /= 2, i--) // convert mapped value into binary string with length "len"
        temp[i] = (tmp%2)?'1':'0';
    temp[len] = '\0';
    strcat(code, temp);
    //cout << temp << " ";
    //cout << strlen(code) << " ";
}

void zigzag(int quant_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE], int zigzag_r[]) // zigzag one block
{
	int pot_x[8*8]={0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,4,5,6,7,7,6,5,4,3,2,1,2,3,4,5,6,7,7,6,5,4,3,4,5,6,7,7,6,5,6,7,7};
	int pot_y[8*8]={0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,7,6,5,4,3,2,3,4,5,6,7,7,6,5,4,5,6,7,7,6,7};
	for (int i=0; i<JPEG_BLOCK_SIZE*JPEG_BLOCK_SIZE; i++)
	{
		zigzag_r[i]=quant_r[pot_x[i]][pot_y[i]];
	}
}

void act_ac_entro(char code[], int zero_num, int val) // entropy code every AC component
{
    char ac_table[16][11][20]={{"1010","00", "01", "100", "1011", "11010", "1111000", "11111000", "1111110110", "1111111110000010", "1111111110000011"},
								{"","1100", "11011", "1111001", "111110110", "11111110110", "1111111110000100", "1111111110000101", "1111111110000110", "1111111110000111", "1111111110001000"},
								{"","11100", "11111001", "1111110111", "111111110100", "1111111110001001", "1111111110001010", "1111111110001011", "1111111110001100", "1111111110001101","1111111110001110"},
								{"","111010", "111110111", "111111110101", "1111111110001111", "1111111110010000", "1111111110010001", "1111111110010010", "1111111110010011", "1111111110010100", "1111111110010101"},
								{"","111011", "1111111000", "1111111110010110", "1111111110010111", "1111111110011000", "1111111110011001", "1111111110011010", "1111111110011011", "1111111110011100", "1111111110011101"},
								{"","1111010", "11111110111", "1111111110011110", "1111111110011111","1111111110100000", "1111111110100001", "1111111110100010", "1111111110100011", "1111111110100100", "1111111110100001"},
								{"","1111011", "111111110110", "1111111110100110", "1111111110100111", "1111111110101000", "1111111110101001", "1111111110101010", "1111111110101011", "1111111110101100", "1111111110101101"},
								{"","11111010", "111111110111", "1111111110101110", "1111111110101111", "1111111110110000", "1111111110110001", "1111111110110010", "1111111110110011", "1111111110110100", "1111111110110101"},
								{"","111111000", "111111111000000", "1111111110110110", "1111111110110111", "1111111110111000", "1111111110111001", "1111111110111010", "1111111110111011", "1111111110111100", "1111111110111101"},
								{"","111111001", "1111111110111110", "1111111110111111", "1111111111000000", "1111111111000001", "1111111111000010", "1111111111000011", "1111111111000100", "1111111111000101", "1111111111000110"},
								{"","111111010", "1111111111000111", "1111111111001000", "1111111111001001", "1111111111001010", "1111111111001011", "1111111111001100", "1111111111001101", "1111111111001110", "1111111111001111"},
								{"","1111111001", "1111111111010000", "1111111111010001", "1111111111010010", "1111111111010011", "1111111111010100", "1111111111010101", "1111111111010110", "1111111111010111", "1111111111011000"},
								{"","1111111010", "1111111111011001", "1111111111011010", "1111111111011011", "1111111111011100", "1111111111011101", "1111111111011110", "1111111111011111", "1111111111100000", "1111111111100001"},
								{"","11111111000", "1111111111100010", "1111111111100011", "1111111111100100", "1111111111100101", "1111111111100110", "1111111111100111", "1111111111101000", "1111111111101001", "1111111111101010"},
								{"","1111111111101011", "1111111111101100", "1111111111101101", "1111111111101110", "1111111111101111", "1111111111110000", "1111111111110001", "1111111111110010", "1111111111110011", "1111111111110100"},
								{"11111111001", "1111111111110101", "1111111111110110", "1111111111110111", "1111111111111000", "1111111111111001", "1111111111111010", "1111111111111011", "1111111111111100", "1111111111111101", "1111111111111110"}};
    char temp[30];
    if ((val == 0)&&(zero_num != 15)) // block end
        strcpy(temp, "1010");
    else
    {
        int len;
        if (val == 0) // zero can't do log2 calculation
            len = 0;
        else
            len = (int)floor(log2(abs(val))) + 1; // length of binary of absolute value
        //strcpy(temp, ac_table[zero_num][len]);
        strcat(code, ac_table[zero_num][len]);
        int unit = (val>0)?val:(val+(2<<(len-1))-1); // mapping (-[2^(len)], -[2^(len-1)]) to (0, 2^(len-1)-1)
        memset(temp, 0, 30*sizeof(char));
        for (int tmp = unit, i = len-1; i>=0; tmp /= 2, i--) // convert mapped value into binary string with length "len"
            temp[i] = (tmp%2)?'1':'0';
        temp[len] = '\0';
    }
    strcat(code, temp);
}

void ac_entro(char code[], int quant_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE]) // entropy code all AC component in one block
{
    int zigzag_r[JPEG_BLOCK_SIZE*JPEG_BLOCK_SIZE];
    int last_no0 = 0;
    int zero_num;
    int val;
    zigzag(quant_r, zigzag_r);
    for (int i=JPEG_BLOCK_SIZE*JPEG_BLOCK_SIZE-1; i>0; i--) // find block end position
        if (zigzag_r[i] != 0)
        {
            last_no0 = i;
            break;
        }
    zero_num = 0;
    for (int i=1; i<=last_no0; i++)
    {
        if ((zigzag_r[i] == 0)&&(zero_num < 15))
            zero_num++;
        else // entropy code when met non-zero element (RLE) or continuous 15 zeros
        {
            val = zigzag_r[i];
            act_ac_entro(code, zero_num, val);
            zero_num = (val==0)?1:0;
        }
    }
    if((last_no0<63)&&(last_no0>=0)) // entropy code block end
        act_ac_entro(code, 0, 0);
}

//int main()
int main(int argc, char *argv[])
{
    FILE *data;
    char ainput[] = "lady.dat"; // default input data file path
    char aoutput[] = "lady_pre.dat"; //default output data file path
    char *pinput = ainput;
    char *poutput = aoutput;
    int width = 256; // default input data width
    int height = 256; // default output data width
    if (argc >= 2)
    {
        pinput = argv[1];
        if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
        {
            cout << "usage: ./jpeg_r [input_file] [output_file] [width] [height]" <<endl;
            exit(0);
        }
    }
    if (argc >= 3)
        poutput = argv[2];
    if (argc >= 5) // resolution of input data
    {
        char *endptr;
        width = strtol(argv[3], &endptr, 10);
        height = strtol(argv[4], &endptr, 10);
    }
    data = fopen(pinput, "r");
    //ifstream input_data (pinput, ios_base::binary);
    if (data == NULL)
    {
        fprintf(stderr, "Input file %s load error!", pinput);
        exit(-1);
    }
    char buffer[width*height]; // store origin data
    fread(buffer, width*height, 1, data); // read the input file
    fclose(data);
    //jpeg_encode(buffer, width, height, poutput);
    char block[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE]; // buffer when dealing with blocks
    double dct_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE]; // store DCT result
    int quant_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE]; // store QUANT result
    char *code;
    code = (char*)malloc(5000000*sizeof(char)); //store entropy coding result
    int PRE_DC = 0;
    //int dc_entro_len = 0;
    int entro_len = 0;
    int total_len = 0;
    memset(code, 0, 5000000*sizeof(char));
    for (int m=0; m<width/JPEG_BLOCK_SIZE; m++)
        for (int n=0; n<height/JPEG_BLOCK_SIZE; n++) // deal with every block
        {
            for (int i=0; i<JPEG_BLOCK_SIZE; i++)
                for (int j=0; j<JPEG_BLOCK_SIZE; j++)
                    block[i][j] = buffer[(m*JPEG_BLOCK_SIZE+i)*width+(n*JPEG_BLOCK_SIZE+j)] - 128;
            dct(block, dct_r);
            quant(dct_r, quant_r);
            int DC = quant_r[0][0]; // DC component
            dc_entro(code, PRE_DC, DC);
            ac_entro(code, quant_r);
            PRE_DC = DC;
        }
    entro_len = strlen(code);
    int count = 0;
    ofstream result(poutput, ios_base::binary);

    for (int i=0; i<entro_len/8; i++) // write every 8 entropy codes as byte into file
    {
        int target = 0;
        for (int j=0; j<8; j++)
        {
            if (code[i*8+j] == '0')
                count = 0;
            else
                count++;
            char tmp = code[i*8+j]-'0';
            target+=tmp<<(7-j);
        }
        char target_r = (char)target;
        result << target_r;
        if(count >= 8) // if continuous 8 '1's, insert a 0 byte into file
        {
            result << '\0';
            count = 0;
        }
    }
    
    int remain;
    if((remain = entro_len % 8) != 0) // pad the last unaligned byte
    {
        int target = 0;
        for (int i=0; i<remain; i++)
        {
            char tmp = code[(entro_len-remain)+i]-'0';
            target+=tmp<<(7-i);
        }
        char target_r = (char)target;
        result << target_r;
    }
    return 0;
}