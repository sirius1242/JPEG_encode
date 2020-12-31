#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>

#define pi 3.1415926
#define JPEG_BLOCK_SIZE 8
using namespace std;
//void dct(char **(origin), double **(dct_r));
//void quant(double **(dct_r), int **(quant_r));
//void jpeg_encode(char buffer[], int width, int height, ofstream result);

void mulmat(double A[], double mat[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE], double B[])
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

void dct(char origin[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE], double dct_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE])
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

void quant(double dct_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE], int quant_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE])
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

int dc_entro(char code[], int PRE_DC, int DC)
{
    char dc_len_table[12][10] = {"00","010", "011","100","101","110","1110","11110","111110","1111110","11111110","111111110"};
    //code[0] = '\0';
    int diff = DC - PRE_DC;
    int leninc = 0;
    int len;
    if (diff == 0)
        len = 0;
    else
        len = (int)floor(log2(abs(diff))) + 1;
    strcat(code, dc_len_table[len]);
    leninc += strlen(dc_len_table[len]);
    int unit = (diff>0)?diff:(diff+(2<<(len-1))-1);
    char temp[16];
    memset(temp, 0, 16*sizeof(char));
    //for (int tmp = unit, i = 0; tmp >= 0; tmp /= 2, i++)
    //    temp[i] = (tmp%2)?'1':'0';
    //for (int tmp = unit, i = len-1; i>=0; tmp /= 2, i--)
    //    temp[i] = (tmp%2)?'1':'0';
    //temp[len] = '\0';
    for (int tmp = unit, i = len-1; i>=0; tmp /= 2, i--)
        temp[i] = (tmp%2)?'1':'0';
    temp[len] = '\0';
    leninc += len;
    strcat(code, temp);
    //cout << temp << " ";
    //cout << strlen(code) << " ";
    return leninc;
}

void zigzag(int quant_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE], int zigzag_r[]){
	int pot_x[8*8]={0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,4,5,6,7,7,6,5,4,3,2,1,2,3,4,5,6,7,7,6,5,4,3,4,5,6,7,7,6,5,6,7,7};
	int pot_y[8*8]={0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,7,6,5,4,3,2,3,4,5,6,7,7,6,5,4,5,6,7,7,6,7};
	for (int i=0; i<JPEG_BLOCK_SIZE*JPEG_BLOCK_SIZE; i++)
	{
		zigzag_r[i]=quant_r[pot_x[i]][pot_y[i]];
	}
}

int act_ac_entro(char code[], int zero_num, int val)
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
    int leninc = 0;
    if ((val == 0)&&(zero_num != 15))
        strcpy(temp, "1010");
    else
    {
        int len;
        if (val == 0)
            len = 0;
        else
            len = (int)floor(log2(abs(val))) + 1;
        //strcpy(temp, ac_table[zero_num][len]);
        strcat(code, ac_table[zero_num][len]);
        leninc += strlen(ac_table[zero_num][len]);
        int unit = (val>0)?val:(val+(2<<(len-1))-1);
        memset(temp, 0, 30*sizeof(char));
        for (int tmp = unit, i = len-1; i>=0; tmp /= 2, i--)
            temp[i] = (tmp%2)?'1':'0';
        temp[len] = '\0';
        leninc += len;
    }
    strcat(code, temp);
    leninc += strlen(temp);
    return leninc;
}

int ac_entro(char code[], int quant_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE])
{
    int zigzag_r[JPEG_BLOCK_SIZE*JPEG_BLOCK_SIZE];
    int last_no0;
    int zero_num;
    int val;
    int leninc = 0;
    zigzag(quant_r, zigzag_r);
    for (int i=JPEG_BLOCK_SIZE*JPEG_BLOCK_SIZE-1; i>0; i--)
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
        else
        {
            val = zigzag_r[i];
            leninc += act_ac_entro(code, zero_num, val);
            zero_num = (val==0)?1:0;
        }
    }
    if((last_no0<63)&&(last_no0>=0))
    {
        leninc += act_ac_entro(code, 0, 0);
    }
    return leninc;
}

//int main()
int main(int argc, char *argv[])
{
    FILE *data;
    char ainput[] = "lady.dat";
    char aoutput[] = "lady_pre.dat";
    char *pinput = ainput;
    char *poutput = aoutput;
    if (argc >= 2)
        pinput = argv[1];
    if (argc >= 3)
        poutput = argv[2];
    data = fopen(pinput, "r");
    //ifstream input_data (pinput, ios_base::binary);
    if (data == NULL)
    {
        fprintf(stderr, "Input file %s load error!", pinput);
        exit(-1);
    }
    int width = 256;
    int height = 256;
    char buffer[width*height];
    fread(buffer, width*height, 1, data);
    fclose(data);
    //jpeg_encode(buffer, width, height, poutput);
    char block[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE];
    double dct_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE];
    int quant_r[JPEG_BLOCK_SIZE][JPEG_BLOCK_SIZE];
    char *code;//[5000000];
    char *result_c;//[5000000];
    code = (char*)malloc(5000000*sizeof(char));
    result_c = (char*)malloc(5000000*sizeof(char));
    int PRE_DC = 0;
    //int dc_entro_len = 0;
    int entro_len = 0;
    int total_len = 0;
    memset(code, 0, 5000000*sizeof(char));
    memset(result_c, 0, 5000000*sizeof(char));
    for (int m=0; m<width/JPEG_BLOCK_SIZE; m++)
        for (int n=0; n<height/JPEG_BLOCK_SIZE; n++) // deal with every block
        {
            for (int i=0; i<JPEG_BLOCK_SIZE; i++)
                for (int j=0; j<JPEG_BLOCK_SIZE; j++)
                    block[i][j] = buffer[(m*JPEG_BLOCK_SIZE+i)*width+(n*JPEG_BLOCK_SIZE+j)] - 128;
            dct(block, dct_r);
            quant(dct_r, quant_r);
            int DC = quant_r[0][0];
            entro_len += dc_entro(code, PRE_DC, DC);
            entro_len += ac_entro(code, quant_r);
            PRE_DC = DC;
        }
    //entro_len = strlen(code);
    cout << entro_len << endl;
    int jpeg_len = 0;
    int count = 0;
    ofstream result(poutput, ios_base::binary);
    //cout << entro_len << endl;
    for (int i=0; i<entro_len; i++)
    {
        result_c[jpeg_len] = code[i];
        jpeg_len++;
        if (code[i] == '0')
            count = 0;
        else
            count++;
        if ((count>=8)&&(jpeg_len%8==0))
        {
            for (int j=0; j<8; j++)
                result_c[jpeg_len++] = '0';
            count = 0;
        }
    }

    for (int i=0; i<jpeg_len/8; i++)
    {
        int target = 0;
        for (int j=0; j<8; j++)
        {
            char tmp = result_c[i*8+j]-'0';
            target+=tmp<<(7-j);
        }
        char target_r = (char)target;
        result << target_r;
    }
    
    if(jpeg_len % 8 != 0)
    {
        int target = 0;
        for (int i=0; i<jpeg_len-(jpeg_len/8*8); i++)
        {
            char tmp = result_c[(jpeg_len/8*8)+i]-'0';
            target+=tmp<<(7-i);
        }
        char target_r = (char)target;
        result << target_r;
    }
    /*
    if (input_data)
    {
        char *buffer = new char[JPEG_BLOCK_SIZE];
        while (input_data.peek() != EOF)
        {
            memset(buffer, 0, sizeof(char)*JPEG_BLOCK_SIZE);
            input_data.read(buffer, )
        }
    }
    */
    return 0;
}