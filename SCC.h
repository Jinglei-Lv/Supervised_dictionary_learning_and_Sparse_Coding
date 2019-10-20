#ifndef SPARSE_COORDINATE_CODING_H
#define SPARSE_COORDINATE_CODING_H

namespace dpl{

static unsigned int myseed;

double getAbs( double value ){
	if( value < 0 )
		return -1*value;
	else 
		return value;
}

double **FeatureInitialization( int featureNumber, int sampleNumber ){

	double **feature = (double**)malloc(sampleNumber*sizeof(double*));
        for( unsigned int i=0; i<sampleNumber; i++ ){
		feature[i] = (double*)malloc(featureNumber*sizeof(double));
		for( unsigned int j=0; j<featureNumber; j++ )
			feature[i][j] = 0;
	}
	return feature;
}

std::vector<int> *NonZeroIndexInitialization( int sampleNumber ){
	std::vector<int> *nonZeroIndex = new std::vector<int> [sampleNumber];
	return nonZeroIndex;
}

double ShrinkageFunction( double value, double theta ){

	if( value < -theta )
		return value+theta;
	else if( value > theta )
		return value-theta;
	else
		return 0;
}

double *Initialize_A_Copy( int featureNumber ){

	double *A_Copy = (double*)malloc(featureNumber*sizeof(double));
	for( unsigned int i=0; i<featureNumber; i++ )
		A_Copy[i]=0;
	return A_Copy;
}

double *Initialize_A( int featureNumber ){

	double *A = (double*)malloc(featureNumber*sizeof(double));
	for( unsigned int i=0; i<featureNumber; i++ )
		A[i]=0;
	return A;
}

void Initialize_A( double *A, double *A_Copy, int featureNumber ){
	for( unsigned int i=0; i<featureNumber; i++ ){
		A[i]=A_Copy[i];
		A_Copy[i]=0;
    	}
}

void Update_A( double *A, double *A_Copy, double *feature, std::vector<int> &nonZeroIndex ){
	for( unsigned int i=0; i<nonZeroIndex.size(); i++ ){
		A[nonZeroIndex[i]] += feature[nonZeroIndex[i]]*feature[nonZeroIndex[i]];
		A_Copy[nonZeroIndex[i]] += feature[nonZeroIndex[i]]*feature[nonZeroIndex[i]];
	}
}

double getNonNegativeFeature( double featureElement, double optimalT ){
	if( featureElement+optimalT>=0 )
		return optimalT;
	else
		return -1*featureElement;
}

int *getRandomIndex( int size ){
 
	std::vector<int> index (size);
	int *data=(int*)malloc(size*sizeof(int));
	for( unsigned int i=0; i<size; i++ )
	        index[i] = i;
 
	for( unsigned int i=0; i<size; i++ ){
    		int randomIndex = rand()%index.size();
        	data[i] = index[randomIndex];
        	index.erase(index.begin()+randomIndex);
    	}	
    	return data;
}

void UpdateFeature( double **Wd, double *sample, double *residuals, double *feature, std::vector<int> &nonZeroIndex, std::vector<int> & nonZeroIndexinitial, double lambda, int layers, int featureNumber, int sampleElementNumber,int fixfeaturebeginindex,int fixfeaturemapnum,double* fixedmap, bool NonNegative ){
    
    	for( unsigned int i = 0; i<sampleElementNumber; i++ ){
  		residuals[i] = -sample[i];
    		for( unsigned int j = 0; j<nonZeroIndex.size(); j++ )
    		        residuals[i] += Wd[i][nonZeroIndex[j]]*feature[nonZeroIndex[j]];
    	}

	nonZeroIndex.resize(0);
	nonZeroIndex=nonZeroIndexinitial;
	int *randomIndex = getRandomIndex(featureNumber );

	for ( unsigned int i = 0; i < featureNumber; i++ ){

        	double optimalT;
        	double derivative = 0;

        	for (unsigned int j = 0;j < sampleElementNumber; j++)
                	derivative += (residuals[j]*Wd[j][randomIndex[i]]);
                
		optimalT = ShrinkageFunction( feature[randomIndex[i]]-derivative, lambda )-feature[randomIndex[i]];

		if( NonNegative ) 
			optimalT = getNonNegativeFeature( feature[randomIndex[i]], optimalT ); 

		feature[randomIndex[i]] += optimalT;

        	if ( optimalT!=0 ){
            		for (unsigned int j = 0;j < sampleElementNumber; j++)
                		residuals[j] += optimalT*Wd[j][randomIndex[i]];
        	}

		if( feature[randomIndex[i]]!=0 && ( randomIndex[i]<fixfeaturebeginindex||randomIndex[i]>=fixfeaturebeginindex+fixfeaturemapnum))
			nonZeroIndex.push_back(randomIndex[i]);
		if(randomIndex[i]>=fixfeaturebeginindex && randomIndex[i]<(fixfeaturebeginindex+fixfeaturemapnum))
		{
			feature[randomIndex[i]]=feature[randomIndex[i]]*fixedmap[randomIndex[i]-fixfeaturebeginindex];
		}

	}	

	for ( unsigned int k = 1; k < layers; k++ ){
		for ( unsigned int i = 0; i < nonZeroIndex.size(); i++ ){
        		double optimalT;
        		double derivative = 0;
        		for (unsigned int j = 0;j < sampleElementNumber; j++)
                		derivative += (residuals[j]*Wd[j][nonZeroIndex[i]]);
                
			optimalT = ShrinkageFunction( feature[nonZeroIndex[i]]-derivative, lambda )-feature[nonZeroIndex[i]];

			feature[nonZeroIndex[i]] += optimalT;

        		if ( optimalT!=0 ){
            			for (unsigned int j = 0;j < sampleElementNumber; j++)
                			residuals[j] += optimalT*Wd[j][nonZeroIndex[i]];
        		}
		}
	}

	nonZeroIndex.resize(0);
	nonZeroIndex=nonZeroIndexinitial;
	for ( unsigned int i = 0; i < featureNumber; i++ ){
		if( feature[i]!=0 && ( i<fixfeaturebeginindex||i>=fixfeaturebeginindex+fixfeaturemapnum))
			nonZeroIndex.push_back(i);
		if( feature[i]==0 && ( i>=fixfeaturebeginindex&&i<fixfeaturebeginindex+fixfeaturemapnum))
		{
			feature[i]=0.0001*fixedmap[i-fixfeaturebeginindex];
		}
	}
	free(randomIndex);
}


void UpdateWd( double **Wd, double *residuals, double *feature, double *A, std::vector<int> &nonZeroIndex, int sampleElementNumber,int fixdicbeginindex,int fixdicnum, bool MinCorrDcDlstatus, float gamma, int featureNumber, double **Wfixdco, double **Wdco ){
		

    	for ( unsigned int i = 0; i < sampleElementNumber; i++ ){
        	for ( unsigned int j = 0; j < nonZeroIndex.size(); j++ )
			{
				if (nonZeroIndex[j]<fixdicbeginindex||nonZeroIndex[j]>=fixdicbeginindex+fixdicnum){
					if(MinCorrDcDlstatus)
					{
						Wdco[i][nonZeroIndex[j]]=0;
						for (unsigned int k=0;k<sampleElementNumber;k++){
							Wdco[i][nonZeroIndex[j]]+=Wfixdco[i][k]*Wd[k][nonZeroIndex[j]];
						}
						Wd[i][nonZeroIndex[j]] = Wd[i][nonZeroIndex[j]] - feature[nonZeroIndex[j]]*residuals[i]*dpl::learningRate(A,nonZeroIndex[j])-gamma*Wdco[i][nonZeroIndex[j]];
					}
					else
					{
						Wd[i][nonZeroIndex[j]] = Wd[i][nonZeroIndex[j]] - feature[nonZeroIndex[j]]*residuals[i]*dpl::learningRate(A,nonZeroIndex[j]);
					}
            		
				}
    	
			}
		}

}

void UpdateWd( double **Wd, double* sample, double *residuals, double *feature, double *A, std::vector<int> &nonZeroIndex, int sampleElementNumber, int fixdicbeginindex, int fixdicnum, bool MinCorrDcDlstatus, float gamma, int featureNumber, double **Wfixdco, double **Wdco ){
		
	for( unsigned int i = 0; i<sampleElementNumber; i++ ){
  			residuals[i] = -sample[i];
    		for( unsigned int j = 0; j<nonZeroIndex.size(); j++ )
    		        residuals[i] += Wd[i][nonZeroIndex[j]]*feature[nonZeroIndex[j]];
    	}

    	for ( unsigned int i = 0; i < sampleElementNumber; i++ ){
        	for ( unsigned int j = 0; j < nonZeroIndex.size(); j++ ){
				    if (nonZeroIndex[j]<fixdicbeginindex||nonZeroIndex[j]>=fixdicbeginindex+fixdicnum)
					{
						if(MinCorrDcDlstatus)
						{
							Wdco[i][nonZeroIndex[j]]=0;
							for (unsigned int k=0;k<sampleElementNumber;k++){
								Wdco[i][nonZeroIndex[j]]+=Wfixdco[i][k]*Wd[k][nonZeroIndex[j]];
							}
							Wd[i][nonZeroIndex[j]] = Wd[i][nonZeroIndex[j]] - feature[nonZeroIndex[j]]*residuals[i]/A[nonZeroIndex[j]]-gamma*Wdco[i][nonZeroIndex[j]];
						}
						else
						{
							Wd[i][nonZeroIndex[j]] = Wd[i][nonZeroIndex[j]] - feature[nonZeroIndex[j]]*residuals[i]/A[nonZeroIndex[j]];
						}
					}
			}
        }

}

void NormalizeWd( double **Wd, std::vector<int> &nonZeroIndex, int sampleElementNumber ){
	for( unsigned int i=0; i<nonZeroIndex.size(); i++ ){
		double sum = 0;
		for( unsigned int j=0; j<sampleElementNumber; j++ )
			sum += Wd[j][nonZeroIndex[i]]*Wd[j][nonZeroIndex[i]];
		sum = sqrt(sum);
		
		if( sum!=0 ){
			for( unsigned int j=0; j<sampleElementNumber; j++ )
				Wd[j][nonZeroIndex[i]] = Wd[j][nonZeroIndex[i]]/sum;
		}
	}
}

void saveFeature( double **feature, char *FeatureFileName, int featureNumber, int sampleNumber ){
	
	printf("Save Features in %s\n", FeatureFileName);

	FILE *fp;
        fp = fopen( FeatureFileName, "w");
        if( fp == NULL ){
		printf("could not find feature file %s\n", FeatureFileName);
            	exit(0);
	}

	for( unsigned int i=0; i<featureNumber; i++ ){
		for( unsigned int j=0; j<sampleNumber; j++)
	        	fprintf(fp, "%.15lf ", feature[j][i]);
		fprintf(fp, "\n");		
	}
	fclose(fp);
}

void saveNonZeroIndex( std::vector<int> *nonZeroIndex, char *IndexFileName, int featureNumber, int sampleNumber ){
	
	printf("Save nonZero index in %s\n", IndexFileName);

	FILE *fp;
        fp = fopen( IndexFileName, "w");
        if( fp == NULL ){
		printf("could not find index file %s\n", IndexFileName);
            	exit(0);
	}

	for( unsigned int i=0; i<sampleNumber; i++ ){
		for( unsigned int j=0; j<nonZeroIndex[i].size(); j++)
	        	fprintf(fp, "%d ", nonZeroIndex[i][j]);
		fprintf(fp, "\n");		
	}
	fclose(fp);
}

void clearFeature( int sampleNumber, double **feature ){
	
	for( unsigned int i=0; i<sampleNumber; i++ )
		free(feature[i]);
	free(feature);
}

double computeLassoResult( double **Wd, double *sample, double *feature,double * resvec, double lambda, int sampleElementNumber, int featureNumber ){
    
	double LassoResult = 0;
	double residuals;
	for( unsigned int i=0; i<sampleElementNumber; i++ ){
		residuals = -sample[i];
		resvec[i]=sample[i];
		for( unsigned int j=0; j<featureNumber; j++ ){
			residuals += Wd[i][j]*feature[j];
		    resvec[i]=resvec[i]-Wd[i][j]*feature[j];
		}
		LassoResult += residuals*residuals;
	}
    
	double sum_feature = 0;
	for( unsigned int j=0; j<featureNumber; j++ )
		sum_feature += getAbs(feature[j]);
    
    	return 0.5*LassoResult+lambda*sum_feature;
}


void calculateError(  double **Wd,  double **sample, double **feature,double ** resvec, double lambda, int sampleNumber, int sampleElementNumber, int featureNumber ) {

	double TotalDecError = 0;
	for( unsigned int t=0; t<sampleNumber; t++ ){
		TotalDecError += computeLassoResult( Wd, sample[t], feature[t],resvec[t], lambda, sampleElementNumber, featureNumber);
	}
	TotalDecError /= sampleNumber;
	std::cout<<"Total Decode Error is "<<TotalDecError<<std::endl;
}

void trainDecoder( double **Wd, double **feature, double **sample,double** fixedmap, double lambda, int layers, int featureNumber, int sampleNumber, int sampleElementNumber, int iterationNumber,int fixdicbeginindex,int fixdicnum, int fixfeaturebeginindex,int fixfeaturemapnum,bool NonNegative,double** residualmat,bool MinCorrDcDlstatus, float gamma, double **Wfixdco, double **Wdco ){

	double *residuals = (double*)malloc(sampleElementNumber*sizeof(double));
	double *A = Initialize_A( featureNumber );
	double *A_Copy = Initialize_A_Copy( featureNumber );
	std::vector<int> *nonZeroIndex = NonZeroIndexInitialization( sampleNumber );
	std::vector<int> *nonZeroIndexinitial = NonZeroIndexInitialization( sampleNumber );
	if (fixfeaturemapnum>0)
	{
		for (int index=0;index<=sampleNumber-1;index++)
		{
			for (int fixit=fixfeaturebeginindex;fixit<=fixfeaturebeginindex+fixfeaturemapnum-1;fixit++ )
			{
				if(feature[index][fixit]!=0)
				{
					nonZeroIndex[index].push_back(fixit);
					nonZeroIndexinitial[index].push_back(fixit);
				}
			}
		}
	}



	srand((unsigned)time(0));
	myseed = (unsigned int) RAND_MAX * rand();	

	std::cout<<"Train decoder"<<std::endl;				
	double ComputionalTime = 0;
	double BeginTime = omp_get_wtime();  
    for( unsigned int it=0; it<iterationNumber; it++ ){
		
		int index = it%sampleNumber;
		if( index==0 )
			Initialize_A( A, A_Copy, featureNumber );


	    Update_A( A, A_Copy, feature[index], nonZeroIndex[index] );
		
		UpdateWd( Wd,sample[index], residuals, feature[index], A, nonZeroIndex[index], sampleElementNumber,fixdicbeginindex,fixdicnum,MinCorrDcDlstatus,gamma,featureNumber, Wfixdco,Wdco);

		NormalizeWd( Wd, nonZeroIndex[index], sampleElementNumber ); 

		UpdateFeature( Wd, sample[index], residuals, feature[index], nonZeroIndex[index],nonZeroIndexinitial[index], lambda, layers, featureNumber, sampleElementNumber,fixfeaturebeginindex,fixfeaturemapnum,fixedmap[index], NonNegative );
		
		Update_A( A, A_Copy, feature[index], nonZeroIndex[index] );
		
		UpdateWd( Wd, residuals, feature[index], A, nonZeroIndex[index], sampleElementNumber, fixdicbeginindex,fixdicnum,MinCorrDcDlstatus,gamma,featureNumber,Wfixdco,Wdco);

		NormalizeWd( Wd, nonZeroIndex[index], sampleElementNumber ); 

		if( it%1==0 )
            		std::cout<<it+1<<" iterations finished"<<std::endl;
	}
	double EndTime = omp_get_wtime();  	
	ComputionalTime += (EndTime-BeginTime);

   	std::cout<<"Finish decoding process:"<<std::endl;
	std::cout<<"Train Decode Time is "<<ComputionalTime<<" seconds."<<std::endl;
	calculateError( Wd, sample, feature,residualmat, lambda, sampleNumber, sampleElementNumber, featureNumber );
	free(A_Copy);
	free(residuals);
	delete [] nonZeroIndex;
}


}

#endif /* Sparse Coordinate Coding */

