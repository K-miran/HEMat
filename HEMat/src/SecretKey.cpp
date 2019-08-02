#include "SecretKey.h"

SecretKey::SecretKey(long logN, long h,  long keydist) {
	long N = 1 << logN;
    
    if(keydist == 0){
        NumUtils::sampleHWT(sx, N, h);
    }
    else{
        NumUtils::sampleZO(sx, N);
    }
	
}
