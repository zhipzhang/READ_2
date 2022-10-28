
/*

    Read CORSIKA output data file (DAT FILE)
    Author zhipzhang
    email  zhipzhang@mail.ustc.edu.cn

*/


#ifdef _THIN_
#define _PARTICLE_LENGTH_ 8
#else
#define _PARTICLE_LENGTH_ 7
#endif

#define _NPARTICLE_ 39
#define _NSUBBLOCK_ 21
#define _SUBBLOCK_LENGTH_ _PARTICLE_LENGTH_*_NPARTICLE_
#define _RECORD_LENGTH_ _SUBBLOCK_LENGTH_ * _NSUBBLOCK_

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <fstream>


void syntax()
{
    printf(" Usage: ./Read infile outfile");
}

int main(int argc, char** argv)
{
    if(argc < 3)
    {
        syntax();
        exit(1);
    }
    FILE* infile = fopen(argv[1], "rb");

    if(!infile)
    {
        std::cout << "Error when opening input file " << argv[1] << std::endl;
        fclose(infile);
    }

    union data
    {
        float data;
        char x[4];
    };
    
    data pdata[_RECORD_LENGTH_];

    TFile* out_file = TFile::Open(argv[2], "recreate");
    int run_num;
    while( !feof(infile))
    {
        fread(&run_num, sizeof(int), 1, infile);
        fread(pdata, 4, 5733, infile);
        fread(&run_num, sizeof(int), 1, infile);
        for(int i = 0; i < 5733 ; i += 273)
        {
            if(strcmp(pdata[i].x , "RUNH") == 0 )
            {
                std::cout << "over";
            }
            else if (strcmp(pdata[i].x, "EVTH") == 0)
            {

            }
            else if (strcmp(pdata[i].x, "EVTE") == 0)
            {

            }
            else if (strcmp(pdata[i].x, "RUNE") == 0)
            {

            }
            else if (strcmp(pdata[i].x, "LONG") == 0)
            {
                std::cout << "yes";
            }
            else
            {
                for( int j = i; j < i + 272; j += 7)
                {
                    if(pdata[j].data = 0)
                    {
                        continue;
                    }
                    
                }
            }

        }
    }
    fclose(infile);



}