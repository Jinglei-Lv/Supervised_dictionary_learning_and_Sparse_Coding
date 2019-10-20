/*
	Sparse Coordinate Coding  version 1.0.2
*/
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <ctime>
#include <iomanip>
#include <string>
#include <omp.h>
#include "DictionaryGeneration.h"
#include "SampleNormalization.h"
#include "LR.h"
#include "SCC.h"

int main(int argc, char* argv[])
{
	if (argc!=21&argc!=19)
	{
		std::cout<<"Parameters: 1 SampleFileName 2 FeatureFileName 3 initializedDictionaryName 4 savedDictionaryName 5 featureNumber 6 sampleElementNumber 7 layers 8 epochNumber 9 lambda 10 DictionaryGenerationState" 
			<<"11 NonNegativeState 12 fixdicbeginindex [0--featureNumber-1] 13 fixdicnum 14 fixfeaturebeginindex [0--featureNumber-1] 15 fixfeaturemapnum 16 fixdicfile 17 fixfeaturefile "
			<<"18 residualfile 19 MinCorrDcDlstatus 20 gamma\n";
		exit(1);
	}
	
	char* SampleFileName =argv[1];// "1.WM.sig.txt";
	char* FeatureFileName = argv[2];//"Feature.txt";
	char* initializedDictionaryName =argv[3]; //"RandomPatchDictionary.txt";
	char* savedDictionaryName =argv[4]; //"Dictionary.txt";
	int featureNumber =atoi(argv[5]); //400;
	int sampleElementNumber =atoi(argv[6]);// 405;
	int layers =atoi(argv[7]);// 3;
	int epochNumber = atoi(argv[8]);//10;
	double lambda=atof(argv[9]);//0.13;
	bool DictionaryGenerationState =atoi(argv[10]);// true;
	bool NonNegativeState = atoi(argv[11]);//true;

	int fixdicbeginindex=atoi(argv[12]);
	int fixdicnum=atoi(argv[13]);
	int fixfeaturebeginindex=atoi(argv[14]);
	int fixfeaturemapnum=atoi(argv[15]);//2;
	char* fixdicfile=argv[16];
	char* fixfeaturefile=argv[17];//"fixfeaturemap.txt";
	char* residualfile=argv[18];

	bool MinCorrDcDlstatus;
	float gamma;
	if(argc==21)
	{
		MinCorrDcDlstatus=atoi(argv[19]);
		gamma=atof(argv[20]);
	}else
	{
		MinCorrDcDlstatus=0;
		gamma=0;
	}


	double **Wd;
	double **fixedD;
	double **feature;
	double **sample;
	double **residualmat;
	double **fixedmap;
	

	int sampleNumber = dpl::getSampleNumber( SampleFileName );
	int iterationNumber = sampleNumber*epochNumber;

	std::cout<<"Number of samples is "<<sampleNumber<<std::endl;
	std::cout<<"Number of samples' element is "<<sampleElementNumber<<std::endl;
	std::cout<<"Number of features is "<<featureNumber<<std::endl;
	std::cout<<"Number of Iterations is "<<iterationNumber<<std::endl;
	std::cout<<"lambda is "<<lambda<<std::endl;
	
	std::cout<<"Begin to read sample."<<std::endl;
	sample = dpl::ReadSample( SampleFileName, sampleNumber, sampleElementNumber );

	dpl::SampleNormalization( sample, sampleNumber, sampleElementNumber );

	std::cout<<"Begin to initialize dictionary."<<std::endl;

	if( DictionaryGenerationState )
		Wd = dpl::GenerateRandomPatchDictionary( featureNumber, sampleElementNumber, sampleNumber, sample );	
	else 
		Wd = dpl::readDictionary( initializedDictionaryName, featureNumber, sampleElementNumber );

	
	if (fixdicnum>0)
	{
		fixedD=dpl::readDictionary(fixdicfile, fixdicnum, sampleElementNumber );
		dpl::DictionaryNormalization(fixdicnum,sampleElementNumber,fixedD);
		for (int i=fixdicbeginindex;i<=fixdicbeginindex+fixdicnum-1;i++)
		{
			for (int j=0;j<=sampleElementNumber-1;j++)
			{
				Wd[j][i]=fixedD[j][i-fixdicbeginindex];
			}
		}
	}


	dpl::DictionaryNormalization( featureNumber, sampleElementNumber, Wd );

	if( DictionaryGenerationState )
		dpl::saveDictionary( featureNumber, sampleElementNumber, Wd, initializedDictionaryName );	
	
	feature = dpl::FeatureInitialization( featureNumber, sampleNumber);
	fixedmap=dpl::ReadSample(fixfeaturefile,sampleNumber,fixfeaturemapnum);

	if (fixfeaturemapnum>0)
	{
		
		for (int i=0;i<=sampleNumber-1;i++)
		{
			for (int j=fixfeaturebeginindex;j<=fixfeaturebeginindex+fixfeaturemapnum-1;j++)
			{
				feature[i][j]=fixedmap[i][j-fixfeaturebeginindex];
			}
		}
	}

	residualmat=dpl::FeatureInitialization( sampleElementNumber, sampleNumber);

	double **Wfixdco;
	double **Wdco;	
//	if (MinCorrDcDlstatus)
//	{
		Wfixdco=dpl::InitializeDictionary(sampleElementNumber,sampleElementNumber);
		
		for (unsigned int i=0;i<sampleElementNumber; i++){
			for ( unsigned int j = 0; j <sampleElementNumber; j++ ){
					Wfixdco[i][j]=0;
					for (unsigned int k=fixdicbeginindex;k<fixdicbeginindex+fixdicnum;k++)
					{
						Wfixdco[i][j]+=Wd[i][k]*Wd[j][k];
					}
			}
		}
		Wdco=dpl::InitializeDictionary( featureNumber, sampleElementNumber );
//	}


	std::cout<<"Begin to train "<<std::endl;
	dpl::trainDecoder( Wd, feature, sample,fixedmap, lambda, layers, featureNumber, sampleNumber,  sampleElementNumber, iterationNumber, fixdicbeginindex, fixdicnum,fixfeaturebeginindex,fixfeaturemapnum, NonNegativeState,residualmat,MinCorrDcDlstatus,gamma, Wfixdco,Wdco);
	std::cout<<"Finish training "<<std::endl;

	dpl::saveDictionary( featureNumber, sampleElementNumber, Wd, savedDictionaryName );	
	dpl::saveFeature( feature, FeatureFileName, featureNumber, sampleNumber );
	dpl::saveFeature( residualmat, residualfile, sampleElementNumber, sampleNumber );
	dpl::clearSample( sampleNumber, sample );
	dpl::clearFeature( sampleNumber, feature );
	dpl::clearDictionary( sampleElementNumber, Wd );
	dpl::clearDictionary( sampleElementNumber,Wfixdco);
	dpl::clearDictionary( sampleElementNumber,Wdco);
	std::cout<<"Hello World!"<<std::endl;
	
	return 0;
}
