#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cmath>
#include <stdint.h>
#include <omp.h>
#include <bitset>


int OLIG_LENGTH = 50;
int seqSimilarity = 75;
int maxseqIdentity;
int maxtilingIdentity;
int mintilingIdentityDiff;
int maxtilingIdentityDiff;
int consecMatch = 15;
double oligc = 0.000001, saltc = 0.075;
int tminterval = 10;
int tmmin=1000, tmmax=-1;
int GC_CONT_MAX= 70;
int GC_CONT_MIN = 30;
int dimerLen = 15, dimerStg = 86;
int hpStem = 6, hpminLoop = 1, hpmaxLoop = 3;
int maxOligs = 1;
int paddingLen = 0;
int tilingLen = OLIG_LENGTH;
int nTiles;
int maxThread=-1;

unsigned long PRIME_MM;
unsigned long PRIME_FH = 3277283;
unsigned long PRIME_IH = 850027;
long gSize;
int POS_SIZE = 8;
int BLOCK = 64;
int nBlocks, nNORs, NORoff, lchars, tNORs, tNORoff;
uint16_t preNOR[65536];
int MAX_NUM_THREAD;
int geneno, olig_total;
long nbytes;
int gseqOff;
bool secStr = false;
bool go = false;
int debugLevel = 0;

const char * seed_Phase1 [] = { //S8W10
    "111010011010111",
    "1100101001000100011011",
    "1111000100010000110101",
    "1100110000101000101101",
    "11100100010010010000111",
    "11010010000000101001000111",
    "1101000100001000001000100111",
    "1101010000010000010000011011"

};

const char * seed_Phase2 [] = { //S8W9
    "111010001101011",
    "1100010010000100010111",
    "1011001010001000010011",
    "1010110000110000001011",
    "11010000001000101000111",
    "110001010000000100101011",
    "110100000100100000010000111",
    "1110010000001000010000010011"
};

//const char * seed_Black [] { //
//    "11111111111",
//    "111111111111111"
//};

using namespace std;

struct OLIGO {
    float Tm;
    int gene;
    int length;
    int uDist;
    int dDist;
    int uPad;
    int dPad;
    unsigned char* sequence;
};

struct entry {
	uint64_t key;
	long *pos;
	int size;
};

unsigned char* substr(unsigned char* dest, unsigned char* src, long srclen, size_t offset, size_t length){
    if(debugLevel > 1) {fprintf(stderr,"Call to substr for %d chars at %d of %ld: %p -> %p\n", length, offset, srclen,src, dest);}
    if(offset < 0 || offset >= srclen){
        return NULL;
    }
    memcpy(dest,&src[offset], length * sizeof(unsigned char));
    dest[length] = '\0';
    return dest;
}

void print_Count(uint8_t *simarray, long* geneLoc){
    if(debugLevel > 0) {fprintf(stderr,"Call to print_Count\n");}
    long oligCount = 0;
    int geneCount = 0;
    for (int i=0; i < geneno; ++i){
        //long oligpergeneCount = 0;
        bool bFound = 0;
        for(long j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j)//(geneLoc[i+1]+geneLoc[i])/2; ++j)
            if (simarray[j] == 0){
                ++oligCount;
                bFound = 1;
            //    ++oligpergeneCount;
            }
        if(bFound) ++geneCount;
        //printf("Gene %d: %d candidates\n", i+1, oligpergeneCount);
    }
    printf("There are %ld candidate oligos from %d genes\n\n",oligCount,geneCount);
}

void print_simArray(uint8_t *simarray){
    if(debugLevel > 0) {fprintf(stderr,"Call to print_simArray\n");}
    for(long i=0; i < gSize; ++i){
        printf("%d",simarray[i]);
        if(i == gSize-1 || i > 0 && i % 50 == 49) printf("\n");
    } 
}

void print_HashTable(entry* T){
    for(int i = 0; i <= PRIME_IH;++i){
        if(T[i].size > 0){
            printf("key)%lu : pos) ",T[i].key);
            for(int j = 0; j < T[i].size; ++j){
                printf("%ld ",T[i].pos[j]);
            }
            printf("\n");
        }
    }
}

/*void print_Blocks(int blocks, int shiftno, uint64_t *hit){
    for(int block = blocks - 1; block >= 0; --block){
        for(int i = BLOCK; i >=0; i-=2){
            while(shiftno-- > 0) {i-=2;}
            
        }

    }
}*/

void print_Time(){
        if(debugLevel > 1) {fprintf(stderr,"Call to print_simArray\n");}
	struct tm *current;
	time_t now;
	
	time(&now);
	current = localtime(&now);
	
	printf("Time is %i:%i:%i\n", current->tm_hour, current->tm_min, current->tm_sec);
}

uint64_t get64bit(long rpos, uint8_t *gCoded){
        if(debugLevel > 5) {fprintf(stderr,"Call to get64bit at %ld\n",rpos);}
	long lpos = rpos - 7;
	if (lpos < 0 && rpos >= 0) 
		lpos = 0;
	uint64_t genomeChunk = 0;
	for (long i = lpos; i <= rpos; ++i) {
    	    genomeChunk = (uint64_t)(genomeChunk << 8 | gCoded[i]);
        }
	return genomeChunk;	
}

uint64_t getvarbit(int byteno, long rpos, uint8_t *Coded){
        if(debugLevel > 4) {fprintf(stderr,"Call to getvarbit at %ld\n",rpos);}
	if (rpos < 0){
		cout << "error!!! in get64bit\n";
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	
	long lpos = rpos - byteno + 1;
	if (lpos < 0 && rpos >= 0)
		lpos = 0;
	uint64_t genomeChunk = 0;
	for (long i = lpos; i <= rpos; ++i)
		genomeChunk = (uint64_t)(genomeChunk << 8 | Coded[i]);
	return genomeChunk;
}

uint64_t get64seed(const char *seedSeq){
        if(debugLevel > 4) {fprintf(stderr,"Call to get64seed\n");}
	int seedWeight=0;
	int seedLength = strlen(seedSeq);
	uint64_t seed=0; 
	for(int j = 0; j < seedLength; ++j) {
		if(seedSeq[j] == '1')
		{
			++seedWeight;
			seed = (uint64_t) (seed << 2 | 3);
		}
		else
			seed = (uint64_t) (seed << 2 | 0);
	}
	return seed;
}

long locate_hashSlot(uint64_t k, entry* T, unsigned long hashSize, bool bNew = false){
    if(debugLevel > 4) {fprintf(stderr,"Call to locate_hashSlot with key %ld\n",k);}
    long i =0, j;
    do{
        j = (k + i) % hashSize;
        if (T[j].key == k || (bNew && T[j].size == 0)) {
            return j;
        }
    } while (++i != hashSize);
    return -1;
}

int search_Hash(uint64_t k, entry *T){
        if(debugLevel > 3) {fprintf(stderr,"Call to search_Hash\n");}
	unsigned long i=0, j;
	do {
		j = (long)((k + i) % PRIME_IH);
		if (T[j].key == k)
			return j;
		++i;
	}while (i!=PRIME_IH && T[j].pos!=NULL);
        if(i == PRIME_IH) printf("Hash table cycled\n");
	return -1;
}

void reverseBlocks(int shiftno, int nblocks, uint64_t* fBlocks, uint64_t* rBlocks){
    if(debugLevel > 4) {fprintf(stderr,"Call to reverseBlocks\n");}
    //Note rBlocks should be initialized to all zeros
    //The rightmost sequence is in the first block, the leftmost in
    //the the last block (ex with 4 char blocks. ACCCGTTA is stored as
    //GTTA ACCC
    //int tmpChars = 0;
    for(int block = 0; block < nblocks; ++block){
        int rblock = nblocks - block - 1;
        uint64_t tmpBlock;
        tmpBlock = ~fBlocks[block];
        for(int i = BLOCK / 2 -1; i >= 0; --i){
            //printf("block %d, rblock %d, i %d, shiftno %d, tmpChars %d\n",block,rblock,i,shiftno,tmpChars);
            rBlocks[rblock] <<= 2;
            rBlocks[rblock] |= tmpBlock & 3;
            tmpBlock >>=2;
        }
    }
    shiftno += 1;
    uint64_t mask = (1 << 2*shiftno) - 1;
    for(int i = 0; i < nblocks; i++){
        rBlocks[i] >>= shiftno * 2;
        if(i < nblocks - 1){
            rBlocks[i] |= (rBlocks[i+1] & mask) << (BLOCK - shiftno * 2);
        }
    }
}

void create_Hash(entry* &hashTable, uint8_t *gCoded, uint64_t seed, int seedlen,long nbytes, unsigned long HashSize){
        if(debugLevel > 1) {fprintf(stderr,"Call to create_Hash\n");}
	//Begin Hashing
	uint64_t window =0, hashKey=0;
	long location = -1;
        int seedBytes = (seedlen % 4) ? seedlen / 4 + 1: seedlen / 4;
        int* lindex = (int*) malloc(HashSize*sizeof(int));
	for (long i=0; i < HashSize; ++i) {
	    hashTable[i].key = 0;
	    hashTable[i].pos = NULL;
	    hashTable[i].size = 0;
            lindex[i] = 0;
	}

        for (int step = 0; step < 2; ++step){
	    window = get64bit(nbytes-1, gCoded);
	    for (long byte = nbytes-1; byte >= seedBytes - 1; --byte) {
	        for (int pass=3; pass >=0 && byte*4 + pass +1 >= seedlen; --pass){
	        	hashKey = (uint64_t) (window & seed);
                        long j = locate_hashSlot(hashKey,hashTable,HashSize,step==0);
                        if(j == -1){
                            printf("Hash cycled in create_Hash\n");
                            exit(EXIT_FAILURE);
                        }
                        if(step == 0){ //Determine how much space is needed
                            if(hashTable[j].size == 0){ //Empty slot
                                hashTable[j].key = hashKey;
                            }
                            hashTable[j].size++;
                        } else { //Assign genome positions
	        	    location = byte*4 + pass; 
                            hashTable[j].pos[lindex[j]++] = location;
                        }
	        	window >>=  2;
	        }
                if(byte >= 8) window |= (((uint64_t) gCoded[byte-8]) << (BLOCK-8));
	    }
            if(step == 0){ //Allocate Space
                for(long i = 0; i < HashSize; ++i) {
                    if(hashTable[i].size > 0){
                        hashTable[i].pos = (long*) malloc(hashTable[i].size*sizeof(long));
                    }
                }
            }
        }
        free(lindex);
}



int identity_Count(uint64_t* ehit, uint64_t* chit){
    if(debugLevel > 3) {fprintf(stderr,"Call to identity_Count\n");}
    int nIdent = 0;
    uint64_t* nolxn;
    nolxn = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
    for(int i = 0; i < nBlocks; ++i){
        nolxn[i] = ~(ehit[i] ^ chit[i]);
    }
    uint16_t* curNor = (uint16_t*) nolxn;
    int i = 0;
    for(i = 0; i < nNORs; ++i){
        nIdent += preNOR[curNor[i]];
    }
    if(NORoff != 0){
        nIdent += preNOR[(uint16_t)(curNor[i] << 2*(8 - NORoff))];
    }
    free(nolxn);
    return nIdent;
}

int getBlocks(uint8_t* gCoded, long pos, uint64_t* hit){
    if(debugLevel > 3) {fprintf(stderr,"Call to getBlocks at %ld\n",pos);}
    int off = pos%4;
    long byte = pos/4;
    //Get Coded sequence
    for(int i = 0; i < nBlocks; ++i){
        hit[i] = get64bit(byte - i * 8, gCoded);
    }
    //Shift off characters until the desired one is rightmost
    for(int i = 0; i < nBlocks; ++i){
        hit[i] >>= 2*(3-off);
        if(i < nBlocks - 1)
            hit[i] |= (hit[i+1] & (uint64_t) ( (1 << 2*(3-off)) - 1)) << (BLOCK - 2*(3-off));
    }
    //If necessary, refill the left bits
    if (OLIG_LENGTH > nBlocks * 32 - (3-off) && byte-nBlocks*8 >=0){
        uint8_t lbyte = gCoded[byte-nBlocks*8];
        hit[nBlocks - 1] |= ((uint64_t)lbyte & (uint64_t) ( (1 << 2*(3-off)) - 1)) 
            << (BLOCK - 2*(3-off));
    }
    //determine how many extra characters are encoded
    return nBlocks * (BLOCK / 2) - OLIG_LENGTH - off;
}



int excessIdentity_Count(uint8_t *gCoded, long pose, long posc,bool reve=false, bool revc=false){
    if(debugLevel > 3) {fprintf(stderr,"Call to excessIdentity_Count at %ld_%d vs %ld_%d\n",pose,reve,posc,revc);}
    uint64_t* hit[2];
    long pos[2] = {pose,posc};
    bool reverse[2] = {reve,revc};
    int kane75;

    for(int i = 0; i < 2; i++){
        hit[i] = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
        int shift = getBlocks(gCoded,pos[i],hit[i]);
        if(reverse[i]){
            uint64_t* tmp;
            tmp = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
            for(int j = 0; j < nBlocks; j++)
                tmp[j] = 0;
            reverseBlocks(lchars,nBlocks,hit[i],tmp);
            for(int j = 0; j < nBlocks; j++)
               hit[i][j] =  tmp[j];
            free(tmp);
        }
    }
    kane75 = identity_Count(hit[0],hit[1]);
    free(hit[0]);
    free(hit[1]);

    return (kane75 - maxseqIdentity);
}

bool checkSubOligos(uint8_t *gCoded, long pose, long posc, uint8_t* simarray, bool reve=false, bool revc=false){
    if(debugLevel > 3) {fprintf(stderr,"Call to checkSubOligo at %ld_%d vs %ld_%d\n",pose,reve,posc,revc);}
    uint64_t* hit[2];
    long pos[2] = {pose,posc};
    bool reverse[2] = {reve,revc};
    int nValid = 0;

    for(int i = 0; i < 2; i++){//0 = e (query) 1 = c (subject)
        hit[i] = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
        int shift = getBlocks(gCoded,pos[i],hit[i]);
        if(reverse[i]){
            uint64_t* tmp;
            tmp = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
            for(int j = 0; j < nBlocks; j++)
                tmp[j] = 0;
            reverseBlocks(lchars,nBlocks,hit[i],tmp);
            for(int j = 0; j < nBlocks; j++)
               hit[i][j] =  tmp[j];
            free(tmp);
        }
    }

    uint64_t* nolxn;
    nolxn = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
    for(int i = 0; i < nBlocks; ++i){
        nolxn[i] = ~(hit[0][i] ^ hit[1][i]);
    }
    int nIdent = 0;
    uint16_t* curNor = (uint16_t*) nolxn;
    int i = 0;
    for(i = 0; i < tNORs; ++i){
        nIdent += preNOR[curNor[i]];
    }
    if(tNORoff != 0){
        nIdent += preNOR[(uint16_t)(curNor[i] << 2*(8 - NORoff))];
    }

    for(int tile = 0; tile < nTiles; ++tile){
       long tilepos[2];

       for(int search = 0; search <2; ++search){
           tilepos[search] = (reverse[search]) ? pos[search] - nTiles + tile + 1 : pos[search] - tile;
       }

       //Eliminate Sub Oligo if excessive identity or previously eliminated
       if(nIdent >= maxtilingIdentity || simarray[tilepos[0]] || simarray[tilepos[1]]){
           for(int search = 0; search < 2; ++search){
               if(reverse[search]){
                   simarray[tilepos[search]] = 1;
               } else {
                   simarray[tilepos[search]] = 1;
               }
           }
       } else{
           nValid++;
       }

       int rBlock = tile / 32;
       int rChar = tile % 32;
       int lBlock = (OLIG_LENGTH - nTiles + tile) / 32;
       int lChar = (OLIG_LENGTH - nTiles + tile) % 32;
       if((nolxn[rBlock] >> 2*(rChar)) & 3 == 3){
           nIdent--;
       }
       if(tile < nTiles - 1 && nolxn[lBlock] >> 2*(lChar + 1) & 3 == 3){
           nIdent++;
       }
    }

    free(nolxn);
    free(hit[0]);
    free(hit[1]);

    return nValid;
}

bool second_Check(uint64_t seed, int seedLen, entry* hashTable, uint8_t* gCoded, long pose, long geneStart,
    long geneEnd, uint8_t* simarray){
    if(debugLevel > 2) {fprintf(stderr,"Call to second_Check at %ld\n",pose);}
    uint64_t hashKey[2]={0,0};
    uint64_t* ehit;
    uint64_t* rhit;
    int w, slot[2], shiftno;
    
    ehit = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
    rhit = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
    shiftno = getBlocks(gCoded,pose,ehit);
    for(int i = 0; i < nBlocks; ++i)
        rhit[i] = 0;
    reverseBlocks(lchars,nBlocks,ehit,rhit);
    int offpos = 0;
    for (int pass=1; pass <= OLIG_LENGTH - seedLen + 1; ++pass){
        hashKey[0] = (uint64_t) (ehit[0] & seed); 
        slot[0] = search_Hash(hashKey[0], hashTable);	
        hashKey[1] = (uint64_t)(rhit[0] & seed);
        slot[1] = search_Hash(hashKey[1], hashTable);
        if(slot[0] == -1 && slot[1] == -1){
            printf("uninitialized hashTable Entry Found, Exiting\n");
            exit(EXIT_FAILURE);
        }
        for(int strand = 0; strand < 2; ++strand){
            if(slot[strand] == -1) continue;
            w=0;
            while (w < hashTable[slot[strand]].size){
                long match_pos = hashTable[slot[strand]].pos[w] + offpos;
            	if (match_pos-offpos > OLIG_LENGTH - 1 && match_pos != pose
                       && (match_pos < geneStart + OLIG_LENGTH - 1 || match_pos > geneEnd)
                       && match_pos < gSize){
                    int off = excessIdentity_Count(gCoded, pose, match_pos,(strand==1));
            	    if (off < mintilingIdentityDiff){ //All suboligos valid	
            	    	++w;
                    } else if(mintilingIdentityDiff < 0 && off < maxtilingIdentityDiff){ //Some suboligos are valid, possibly all
                       if(checkSubOligos(gCoded,pose,match_pos,simarray,(strand==1))){
                            ++w;
                        } else{
                            free(ehit);
                            free(rhit);
                            return false;
                        }
                    } else { //No suboligos are valid
                        w = hashTable[slot[strand]].size;
                        for(long i = off; i >= 0; --i){
                            //printf("mp: %ld, i: %ld, off: %ld, pose: %ld\n",match_pos, i, off, pose);
                            if(match_pos - i > OLIG_LENGTH -1) simarray[match_pos - i] = 1;
                            if(match_pos + i < gSize) simarray[match_pos + i] = 1;
                            if(pose - i > OLIG_LENGTH - 1) simarray[pose - i] = 1;
                            if(pose + i < gSize) simarray[pose + i] = 1;
                        }

                        free(ehit);
                        free(rhit);
            	        return  false;
                    }
            	}
            	else{
                    ++w;
                }
            }
        }
        for(int i = 0; i < nBlocks; ++i){
            ehit[i] >>= 2;
            rhit[i] >>= 2;
            if(i < nBlocks - 1) {
                ehit[i] |= (ehit[i+1] & (uint64_t)  (3)) << (BLOCK - 2);
                rhit[i] |= (rhit[i+1] & (uint64_t)  (3)) << (BLOCK - 2);
            }
        }
        ++offpos;
    }
    free(ehit);
    free(rhit);
    return true;
}

uint64_t* encodeOLIGO(uint64_t* Coded,OLIGO x){
    int nblocks = x.length/(BLOCK / 2);
    if(x.length % (BLOCK / 2)) nblocks++;
    if(!Coded) return 0;
    for(int block = nblocks -1; block >= 0; --block){
        uint64_t temp = 0;
        int stpoint = (block == 0) ? x.uPad : x.uPad + (BLOCK / 2) * (block - nblocks) + x.length;
        int edpoint = x.uPad + (BLOCK / 2) * (1+ block - nblocks) + x.length - 1;
        //printf("s: %d\te: %d\n",stpoint,edpoint);
        for(int i = stpoint; i <= edpoint; i++) {
            switch (x.sequence[i]) {
                case 'A':
                    temp = (uint64_t) (temp << 2 | 0);
                    break;
                case 'C':
                    temp = (uint64_t) (temp << 2 | 1);
                    break;
                case 'G':
                    temp = (uint64_t) (temp << 2 | 2);
                    break;
                case 'T':
                    temp = (uint64_t) (temp << 2 | 3);
                    break;
                default:	
                    printf("Error! %c in oligo at %d is not a valid nucleotide: Encoded as A\n",x.sequence[i] , x.uDist);
                    temp = (uint64_t) (temp << 2 | 0);					
            }
        }
        Coded[nblocks - block - 1] = temp;
    }
    return Coded;
}

bool blacklist_Check(uint64_t seed, int seedLen, entry* hashTable, long bLen, uint8_t* bCoded, OLIGO olig){
    if(debugLevel > 2) {fprintf(stderr,"Call to blacklist_Check of oligo at %d in %d\n",olig.uDist,olig.gene);}
    uint64_t hashKey[2]={0,0};
    uint64_t* ehit;
    uint64_t* rhit;
    uint64_t window[2];
    int w, slot[2], shiftno;
    int threshold = olig.length * seqSimilarity / 100;
    

    ehit = (uint64_t *) malloc (nBlocks * sizeof(uint64_t));
    if(!encodeOLIGO(ehit,olig)){
        cout << "Could not encode OLIGO\n";
        exit(EXIT_FAILURE);
    }

    rhit = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
    for(int i = 0; i < nBlocks; ++i)
        rhit[i] = 0;
    reverseBlocks(lchars,nBlocks,ehit,rhit);
    int offpos = 0;
    window[0] = ehit[0];
    window[1] = rhit[0];

    for (int pass=1; pass <= olig.length - seedLen + 1; ++pass){
        for(int strand = 0; strand < 2; strand++){
            hashKey[strand] = (uint64_t) (window[strand] & seed); 
            slot[strand] = search_Hash(hashKey[strand], hashTable);	
        }
        if(slot[0] == -1 && slot[1] == -1){
            continue;
        }

        for(int strand = 0; strand < 2; ++strand){
            if(slot[strand] == -1) continue;
            w=0;
            while (w < hashTable[slot[strand]].size){
                long match_pos = hashTable[slot[strand]].pos[w] + offpos;
            	if (match_pos-offpos > olig.length - 1 && match_pos < bLen){
                    uint64_t* bhit;
                    bhit = (uint64_t*) malloc(nBlocks * sizeof(uint64_t));
                    getBlocks(bCoded,match_pos,bhit);
                    int count = identity_Count((strand == 0) ? ehit : rhit,bhit);
                    free(bhit);
            	    if (count < threshold && count < consecMatch){ //No Match
            	    	++w;
                    } else if(count >= threshold){ //Similarity Match
                        free(ehit);
                        free(rhit);
                        return false;
                    } else { //Potential Consecutive Identity Match
                        //TODO: Implement a check here
                        ++w;
                    }
                } else {
                    ++w;
                }
            }
        }
        for(int strand = 0; strand < 2; strand++){
            window[strand] >>= 2;
            int curBlock = (pass + seedLen - 1) / (BLOCK / 2);
            int posinBlock = (pass + seedLen - 1) % (BLOCK / 2);
            int adjust = BLOCK - 2 *(posinBlock +1);
            if(curBlock < nBlocks) {
                window[strand] |= (((strand == 0) ? ehit : rhit)[curBlock] & (uint64_t) (3 << (posinBlock * 2))) << adjust;
            }
        }
        ++offpos;
    }

    free(ehit);
    free(rhit);
    return true;
}

int kane_75(uint8_t *gCoded, long pose,  long posc, uint8_t *simarray, int curseqc, long* geneLoc){
    if(debugLevel > 3) {fprintf(stderr,"Call to kane_75 at %ld vs %ld\n",pose,posc);}
    int elims = 0;
    long starte = pose, startc = posc;
    int off = -1, firstoff = -1;
    int dir = -1;
    int curseqe =curseqc;
    while(pose < geneLoc[curseqe+1])
        ++curseqe;
    int startseqc = curseqc, startseqe = curseqe;
    do{
        //Calculate Excessive similarity
        off = excessIdentity_Count(gCoded, pose, posc);
        //Eliminate Position and positions on either
        //side for each excess identity
        if(firstoff == -1) firstoff = off;
        for(int i = off; i >= 0; --i){
            if(posc - i > OLIG_LENGTH - 2) simarray[posc-i] = 1;
            simarray[posc+i] = 1;
            simarray[pose-i] = 1;
            if(pose + i > gSize - 1) simarray[pose+i] = 1;
            elims += 4;
        }
        //Slide to the next position
        if(off >= 0){
            posc += dir * (off+1);
            pose += dir * (off+1);
            if(posc < geneLoc[curseqc]) {--curseqc;}
            if(posc >= geneLoc[curseqc+1]) {++curseqc;}
            if(pose < geneLoc[curseqe]) {--curseqe;}
            if(pose >= geneLoc[curseqe]) {++curseqe;}
            //If the left is hit or the genes match, start sliding right
            if(posc < OLIG_LENGTH - 1 || (curseqc == curseqe && dir == -1)){
                dir = 1;
                posc = startc + (firstoff+1);
                pose = starte + (firstoff+1);
                curseqe = startseqe;
                curseqc = startseqc;
                if(posc >= geneLoc[curseqc+1]) {++curseqc;}
                if(pose >= geneLoc[curseqe]) {++curseqe;}
            }
        }
    }while(off >= 0 && pose < gSize -1 && curseqe != curseqc);
    return elims;
}

void kane_15(long pose,  long posc, uint8_t *simarray,long assoc, long* assocList){
    if(debugLevel > 3) {fprintf(stderr,"Call to kane_15 at %ld vs %ld\n",pose,posc);}
    for(int j = -2; j < assoc; ++j){
        long pos;
        if(j == -2) pos = posc;
        if(j == -1) pos = pose;
        else if(j >= 0) pos = assocList[j]; 
        long ranger = pos + tilingLen - consecMatch;
        ranger =  ranger < gSize ? ranger : gSize -1;
        for(long i = pos; i <= ranger; ++i){
            simarray[i] = 1;
        }
    }
    //fprintf(stderr,"olig at %ld matches at %ld\n",pose,posc);
}

void gc_Content(uint8_t* gCoded, uint8_t *simarray,long nbytes){
    if(debugLevel > 1) {fprintf(stderr,"Call to gc_Content\n");}
    long gccount = 0;
    int loff = 4;
    int roff = 4;
    long pos = gSize;
    int lbyte = nbytes - 1;
    int rbyte = nbytes - 1;
    uint8_t lseq = gCoded[lbyte];
    uint8_t rseq = gCoded[rbyte];
    long counted = 0;
    //construct of scan tilingLen nucleotides, then slide this window along the genome
    //'Window' is composed of a left buffer (the nucleotides about to
    //enter the window, and a right buffer (those about to leave)
    while(lbyte >= 0){ //stop when the window would pass outside the genome
        if((lseq & 3) == 1 || (lseq & 3) == 2) ++gccount; //increment count if G or C enters window on the left
        // shift next nucleotide to right of buffer, refill as necessary
        lseq >>= 2; 
        if(--loff == 0) {
            if(lbyte-- > 0){  
                lseq=gCoded[lbyte];
                loff = 4;
            }
        }
        //Once an oligos worth of nucleotides have been scanned
        if(++counted >= tilingLen){
            --pos;
            //Eliminate positions with GC's out of range
            int temp1 = tilingLen * GC_CONT_MIN / 100.0;
            int temp2 = tilingLen * GC_CONT_MAX / 100.0 + 1;
            if(gccount < temp1 || gccount >= temp2)
                simarray[pos] = 1;
            // decrement count if G or C leaves window on the right
            if((rseq & 3) == 1 || (rseq & 3) == 2) --gccount;
            rseq >>= 2;
            if(--roff == 0) {
                rseq=gCoded[--rbyte];
                roff = 4;
            }
        }
    }
}

void melting_Temp(uint8_t *gCoded, uint8_t *simarray, long* geneLoc, float *tmarray){
        if(debugLevel > 1) {fprintf(stderr,"Call to melting_Temp\n");}
	long byte;
        int off;
	uint64_t halfl2=0, halfl=0, halfr=0, nolxnr=0, nolxnl=0, nolxnl2=0; 		
	const float enthalpy[] = {-7.9,-8.4,-7.8,-7.2,-8.5,-8.0,-10.6,-7.8,-8.2,-9.8,-8.0,-8.4,-7.2,-8.2,-8.5,-7.9};
	const float  entropy[] = {-22.2,-22.4,-21.0,-20.4,-22.7,-19.9,-27.2,-21.0,-22.2,-24.4,-19.9,-22.4,-21.3,-22.2,-22.7,-22.2};
	float enth, entr, Tm;	
	int length = tilingLen;
        
        //printf("melting_temp initialized\n");
	
	for (int qos=0; qos < geneno; ++qos){
	    for(long pos = geneLoc[qos] + length - 1; pos < geneLoc[qos+1]; ++pos){
		if (simarray[pos] == 0){
		    byte = pos/4; off = pos%4; 
		    halfr = get64bit(byte, gCoded); halfl = get64bit(byte - 8, gCoded); 
		    if (byte - 16 >= 0)
		    	halfl2 = get64bit(byte - 16, gCoded);
		    
		    halfr >>= 2*(3-off);
		    halfr |= (halfl & (uint64_t) ( (1 << 2*(3-off)) - 1)) << (BLOCK - 2*(3-off));			
		    halfl >>= 2*(3-off);	
		    halfl |= (halfl2 & (uint64_t) ( (1 << 2*(3-off)) - 1)) << (BLOCK - 2*(3-off));			
		    halfl2 >>= 2*(3-off);
		    
		    ////////////////////////////////// BEGIN: Computing Tm		
		    enth = 0; entr = 0;	
		    nolxnr = halfr;
		    nolxnl = halfl;
		    nolxnl2 = halfl2;
		    
		    if ((nolxnr & 3) == 0 || (nolxnr & 3) == 3){enth+=2.3; entr+=4.1;}
		    if ((nolxnr & 3) == 1 || (nolxnr & 3) == 2){enth+=0.1; entr+=-2.8;}
		    
		    for (int i=1; i<=32; ++i) {
		    	if (i !=32){
		    		enth+=enthalpy[nolxnr & 15];
		    		entr+=entropy[nolxnr & 15];
		    	}
		    	else{
		    		enth+=enthalpy[((nolxnl & 3) << 2) | (nolxnr & 3)];
		    		entr+=entropy[((nolxnl & 3) << 2) | (nolxnr & 3)];
		    	}					
		    	nolxnr >>= 2;
		    }	
		    int up1 = (length > 64? 64:length);
		    for (int i=33; i <= up1; ++i) {
		    	if (i  != up1 ){
		    		enth+=enthalpy[nolxnl & 15];
		    		entr+=entropy[nolxnl & 15];							
		    	}
		    	else{				
		    		if (length == up1){
		    			if ((nolxnl & 3) == 0 || (nolxnl & 3) == 3){enth+=2.3; entr+=4.1;}
		    			if ((nolxnl & 3) == 1 || (nolxnl & 3) == 2){enth+=0.1; entr+=-2.8;}
		    		}
		    		
		    		if (length > up1){
		    			enth+=enthalpy[((nolxnl2 & 3) << 2) | (nolxnl & 3)];
		    			entr+=entropy[((nolxnl2 & 3) << 2) | (nolxnl & 3)];										
		    		}
		    	}
		    	nolxnl >>= 2;
		    }	
		    for (int i=65; i <= length; ++i) {		
		    	if (i != length){
		    		enth+=enthalpy[nolxnl2 & 15];
		    		entr+=entropy[nolxnl2 & 15];
		    	}		
		    	else{
		    		if ((nolxnl2 & 3) == 0 || (nolxnl2 & 3) == 3){enth+=2.3; entr+=4.1;}
		    		if ((nolxnl2 & 3) == 1 || (nolxnl2 & 3) == 2){enth+=0.1; entr+=-2.8;}
		    	}
		    	
		    	nolxnl2 >>= 2;
		    }
		    
		    
		    
		    /*enth*=1000;
		     Tm = enth / (entr + 1.987*log(oligc/4)) + 12.0 * log10(saltc) -273.15;
		     */
		    
		    enth=enth*1000;
		    Tm = enth / (entr + 1.987*log(oligc/4)) + 12.0 * log10(saltc) -273.15;
		    
		    
		    ////////////////////////////////// END: Computing Tm
		    tmarray[pos] = Tm;
		}
	    }
	}
        //printf("\npost- for loop\n");
	
	if (tmmax == -1 && tmmin == 1000){
		int count = 0;	double sumTm = 0.0, maxTm=0.0, minTm=100.0;
		for (long i=0; i < geneno; ++i)
			for(long j = geneLoc[i] + length - 1; j < geneLoc[i+1]; ++j)
				if (simarray[j] == 0){
					++count;
					sumTm += tmarray[j]; 
					if (tmarray[j]>maxTm) maxTm=tmarray[j];
					if (tmarray[j]<minTm) minTm=tmarray[j];
				}		
                if(count > 0){
		    float avtm = sumTm / count;
		    printf("max_Tm: %.3f min_Tm: %.3f avg_Tm: %.3f\n",maxTm,minTm,avtm);
		    
		    int uppTm,lowTm; uppTm=(int) (100*maxTm); lowTm=(int) (100*minTm); 
		    //printf("uppTm: %d \t lowTm: %.d \n",uppTm,lowTm);	
		    // Computing Upper Bound Considering the Tm Range 
		    bool **tmFreq = (bool **) malloc(geneno * sizeof(bool *));
		    for (int i=0; i < geneno; ++i)
		    	tmFreq[i] = (bool *) malloc ((uppTm+1)*sizeof(bool));
		    int *sumPm = (int *) malloc (geneno * sizeof(int));
		    
		    for (int i=0; i < geneno; ++i){
		    	for (int j=0; j <= uppTm; ++j)
		    		tmFreq[i][j]=false;
		    	sumPm[i] = 0;
		    }
		    
		    for (int i=0; i < geneno; ++i)
		    	for(long j = geneLoc[i] + length - 1; j < geneLoc[i+1]; ++j)
		    		if (simarray[j] == 0)
		    			tmFreq[i][(int)(100*tmarray[j])]=true;
		    
		    int tmRange = 100*tminterval, upperBound, maxUbound = 0, maxLindex = -1, maxRindex = -1;
		    
		    upperBound = 0;
		    for (int k=0; k < geneno; ++k){
		    	for (int j=uppTm; j > uppTm - tmRange; --j)		
		    		if (tmFreq[k][j]) ++sumPm[k];
		    	if(sumPm[k] > 0) 
		    		++upperBound;
		    }
		    //printf("Upper Bound in the range [%d,%d) is: %d\n", uppTm - tmRange+1, uppTm+1, upperBound);
		    if (upperBound > maxUbound) {maxUbound = upperBound; maxLindex = uppTm - tmRange+1; maxRindex = uppTm;}
		    
		    for (int i=uppTm-1; i >= lowTm + tmRange - 1; --i){
		    	for (int k=0; k < geneno; ++k){
		    		int deltaPm = 0;
		    		if(tmFreq[k][i+1]) --deltaPm;
		    		if(tmFreq[k][i-tmRange+1]) ++deltaPm;
		    		if(deltaPm==1){
		    			if(sumPm[k]==0) ++upperBound;
		    		}
		    		else if (deltaPm==-1){
		    			if(sumPm[k]==1)	--upperBound;
		    		}
		    		sumPm[k]+=deltaPm;
		    	}
		    	//printf("Upper Bound in the range [%d,%d) is: %d\n", i - tmRange+1, i+1, upperBound);
		    	if (upperBound > maxUbound) {maxUbound = upperBound; maxLindex = i - tmRange+1; maxRindex = i;}
		    }
		    printf("Optimal Tm range of length %d is [%.2f, %.2f)\n", tminterval, maxLindex/100.0, (maxRindex+1)/100.0);
		    
		    for (int i=0; i < geneno; ++i)
		    	for(long j = geneLoc[i] + length - 1; j < geneLoc[i+1]; ++j)
		    		if (simarray[j] == 0) 
		    			if(tmarray[j] >=  1.0*maxRindex/100.0+0.01  || tmarray[j] < 1.0*maxLindex/100.0 )
		    				simarray[j] = 1;
		    
		    for(int i=0 ; i < geneno ; ++i)
		    	free(tmFreq[i]);//delete tmFreq[i];
		    free(tmFreq);//delete [] tmFreq;
		    free(sumPm);//delete [] sumPm;
                }
                else{
                    printf("No oligos found to calculate optimal Tm\n");
                }
	}
	else{
            double maxTm = -1, minTm = 1000;
	    for (long i = length - 1; i < gSize; ++i){
		if ((simarray[i] == 0) &&
                    ((tmmax != -1 && tmarray[i] > tmmax) || (tmmin != 1000 && tmarray[i] < tmmin))){
		    simarray[i] = 1;
                } else if(simarray[i] == 0){
                    if(tmarray[i] > maxTm) maxTm = tmarray[i];
                    if(tmarray[i] < minTm) minTm = tmarray[i];
                }
            }
            printf("Requested Tm Interval: ");
            if(tmmax != -1 && tmmin != 1000){
                printf("[%d, %d]",tmmin,tmmax);
            } else if(tmmax == -1){
                printf("[%d, +inf)",tmmin);
            } else {
                printf("(-inf, %d]",tmmax);
            }
            printf("\nObserved Tm Interval: [%.2f,%.2f]\n",minTm, maxTm);
	}
        //printf("End of melting temp function\n");
}

long hash_Search(uint64_t k, uint64_t *T, uint64_t s, unsigned long PRIME){
        if(debugLevel > 4) {fprintf(stderr,"Call to hash_Search\n");}
	long i=0, j;
	do {
		j = (k + i) % PRIME;
		if (T[j] == k)
			return j;
		++i;
	}while (i!=PRIME && T[j]!=s+1);
        if(i==PRIME) printf("Hash Table Cycled\n");
	return -1;
}

long hash_Insert(uint64_t k, uint64_t *T, uint64_t s, unsigned long PRIME){
        if(debugLevel > 4) {fprintf(stderr,"Call to hash_Insert\n");}
    	long i=0, j;
	do {
		j = (k + i) % PRIME;
		if (T[j] == s+1){
			T[j] = k;
			return j;
		}
		else
			++i;
	}while (i!=PRIME);
        if(i == PRIME) printf("Hash Table Cycled\n");
	return 0;
}			


unsigned long identity_Check(uint8_t *gCoded, uint8_t *simarray, long *geneLoc, long nbytes){
        if(debugLevel > 0) {fprintf(stderr,"Call to identity_Check\n");}
	//Begin Coding "111111111111111"
	uint64_t seed=0; 
	for(int j = 0; j < consecMatch; ++j) 
		seed = (uint64_t) (seed << 2 | 3);
	
        //printf("identity_Check: Begun Coding\n");
	//Begin Hashing
	uint64_t window = 0, hashKey[2]= {0,0};
	uint8_t nexwin=0;
	uint64_t *hashTable = (uint64_t *) malloc(PRIME_MM*sizeof(uint64_t));
	long *posTable = (long *) malloc (PRIME_MM*sizeof(long));
        long *assocCount = (long *) malloc (gSize*sizeof(long));
        long **assocTable = (long **) malloc (gSize*sizeof(long*));
	for (long i=0; i < PRIME_MM; ++i) {
		hashTable[i] = seed + 1;
		posTable[i] = -1;
	}
        for (long i =0; i < gSize; ++i){
            assocCount[i] = 0;
            assocTable[i] = NULL;
        }

        long maxassociations = 0, associations = 0;
        //printf("identity_check: Begun Hashing\n");
	
	int SEED_BYTE = consecMatch % 4 == 0 ?  consecMatch / 4 : consecMatch / 4  + 1;
	long rpos = nbytes-1;
	window = getvarbit(SEED_BYTE, rpos, gCoded);
	rpos -= SEED_BYTE;
	int nexwinSize =   SEED_BYTE * 4 - consecMatch;
	if (nexwinSize == 0){
		nexwin = gCoded[rpos];
		--rpos;
		nexwinSize = 4;
	}
	else{
	    nexwin = (window >> 2*consecMatch) & ((1 << 2*nexwinSize)-1); //         
        }
        //printf("indentity_Check: Windows!\n");        
	hashKey[0] = (uint64_t) (window & seed);
        reverseBlocks(BLOCK/2 - consecMatch,1,hashKey+0,hashKey+1);
        hashKey[1] &= seed;
	int curSeq = geneno - 1;
	for (long location = gSize-1; location >= consecMatch-1; --location) {
                if(location < geneLoc[curSeq]){--curSeq;}
                long searchRes[2];
                searchRes[0] = hash_Search(hashKey[0], hashTable, seed, PRIME_MM);
                searchRes[1] = hash_Search(hashKey[1], hashTable, seed, PRIME_MM);
                if(searchRes[0] != -1 || searchRes[1] != -1){
                    for(int strand = 0; strand < 2; ++strand){
		        if (searchRes[strand] != -1){ 
                            long pose = posTable[searchRes[strand]];
                            if(pose >= geneLoc[curSeq] && pose < geneLoc[curSeq + 1]){//match within curSeq
                                if(++assocCount[pose] % consecMatch == 1){
                                    int newsize = (assocCount[pose]/consecMatch +1)*consecMatch;
                                    assocTable[pose] = (long*) realloc(assocTable[pose],newsize*sizeof(long));
                                    for(int i = assocCount[pose]; i < newsize;++i){assocTable[pose][i] = -1;}
                                }
                                assocTable[pose][assocCount[pose]-1] = location;
                                if(++associations > maxassociations) maxassociations = associations;
                            } else{
		                kane_15(pose, location, simarray, assocCount[pose],assocTable[pose]);
                                associations -= assocCount[pose];
                                assocCount[pose] = 0;
                                free(assocTable[pose]);
                                assocTable[pose] = NULL;
                            }
                        }
                    }
                } else{ //If neither the forward or reverse key are in the hash table, add the forward key
                    posTable[hash_Insert(hashKey[0], hashTable, seed, PRIME_MM)] = location;
                }
                // Shift the rightmost character off and fill in the
                // left bits that are now coming into the matching area
		hashKey[0] = (hashKey[0] >> 2) | ( ((uint64_t)(nexwin & 3)) << (consecMatch*2 - 2) );
                // Same operation in reverse with the reverse complement
                hashKey[1] = (hashKey[1] << 2) | ((uint64_t)(~nexwin & 3));
                hashKey[1] &= seed;
		nexwin >>= 2;  --nexwinSize;
		if (nexwinSize == 0 && rpos >= 0){
			nexwin = gCoded[rpos];
			--rpos;
			nexwinSize = 4;
		}
	}
	
	free (hashTable);//delete [] hashTable;
	free(posTable);//delete [] posTable;
        free(assocCount);
        for(long i=0; i < gSize;++i){
            free(assocTable[i]);
        }
        free(assocTable);
        return maxassociations;
}

bool secondary_Structure(uint8_t *gCoded, long pose){
        if(debugLevel > 1) {fprintf(stderr,"Call to secondary_Structure\n");}
	long ebyte = pose/4;
        int eoff = pose%4;
	uint64_t ehitl, ehitm, ehitr, dcwindl=0, rcwindl=0, dcwindm=0, rcwindm=0, dcwindr=0, rcwindr=0, tpwind=0;	
	if(OLIG_LENGTH > 64){
		ehitr = get64bit(ebyte, gCoded); ehitm = get64bit(ebyte - 8, gCoded); ehitl = get64bit(ebyte - 16, gCoded);
		ehitr >>= 2*(3-eoff);
		ehitr |= (ehitm & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitm >>= 2*(3-eoff);
		ehitm |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitl >>= 2*(3-eoff);
		
		dcwindr = ehitr;
		tpwind = ~dcwindr;
		for(int i=1; i <= 32; ++i){
			rcwindl = (rcwindl << 2) | (tpwind & 3);
			tpwind >>= 2;	
		}
		
		tpwind = 0;
		dcwindm = ehitm;
		tpwind = ~dcwindm;
		for(int i=33; i <= 64; ++i){
			rcwindm = (rcwindm << 2) | (tpwind & 3);
			tpwind >>= 2;	
		}
		
		tpwind = 0;
		dcwindl = ehitl;
		tpwind = ~dcwindl;
		for(int i=65; i <= OLIG_LENGTH; ++i){
			rcwindr = (rcwindr << 2) | (tpwind & 3);
			tpwind >>= 2;
		}
		
		rcwindr |= ((uint64_t)(rcwindm & (((uint64_t)1 << 2*(32-(OLIG_LENGTH % 32)))-1)) << 2*(OLIG_LENGTH % 32));
		rcwindm >>= 2*(32-(OLIG_LENGTH % 32));
		
		rcwindm |= ((uint64_t)(rcwindl & (((uint64_t)1 << 2*(32-(OLIG_LENGTH % 32)))-1)) << 2*(OLIG_LENGTH % 32));
		rcwindl >>= 2*(32-(OLIG_LENGTH % 32));
		
		dcwindl &= (( (uint64_t)1 << 2*(OLIG_LENGTH % 32))-1);
		/*printf("dcwind: %16llX \t %16llX \t %16llX \n", dcwindl, dcwindm, dcwindr);
		 printf("rcwind: %16llX \t %16llX \t %16llX \n", rcwindl, rcwindm, rcwindr);*/
		
		// start checking for DIMER
		uint64_t dimerWin=0, dimerTmp=0, dimerSeed=0, dimerKey=0;
		dimerWin = dcwindr;
		dimerTmp = dcwindm;
		dimerSeed = ((uint64_t)1 << 2*dimerLen) - 1;
		//printf("dimerWin: %16llX \n", dimerWin);
		uint64_t *dimerHash = (uint64_t *) malloc((OLIG_LENGTH-dimerLen+1)*sizeof(uint64_t));
		int dtmpCount=0;
		for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i){
			dimerKey = dimerWin & dimerSeed;
			dimerHash[i] = dimerKey;
			dimerWin = (dimerWin >> 2) | ( (uint64_t)(dimerTmp & 3) << (BLOCK-2)); 
			dimerTmp >>= 2; ++dtmpCount;
			if(dtmpCount==32) dimerTmp = dcwindl;
		}
		/*for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i)
		 printf("dimerHash[%d]: %16llX \n",i, dimerHash[i]);*/
		dimerWin = rcwindr;
		dimerTmp = rcwindm;
		for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i){
			dimerKey = dimerWin & dimerSeed;
			//printf("dimerKey: %16llX \n", dimerKey);
			// start core
			for (int j=0; j <= OLIG_LENGTH-dimerLen; ++j){
				int strGncy = 0;
				uint64_t xnorWins = ~(dimerKey ^ dimerHash[j]);
				xnorWins &= dimerSeed;
				for (int k=0; k < (dimerLen-1)/8+1; ++k)
					strGncy += preNOR[((uint16_t) (xnorWins >> 16*k))];
				//cout << (float)strGncy / (float)(dimerLen) << "\n";
				if ( (float)strGncy / (float)(dimerLen) > (dimerStg*1.0)/100.0) {
					/*printf("secondary structure found!\n");
					 printf("dcwind: %16llX \t %16llX \t %16llX %16llX \n", dcwindl, dcwindm, dcwindr, dimerHash[j]);
					 printf("rcwind: %16llX \t %16llX \t %16llX %16llX \n", rcwindl, rcwindm, rcwindr, dimerKey);*/
					free(dimerHash);//delete [] dimerHash;
					return false;
				}
			}
			// end core
			dimerWin = (dimerWin >> 2) | ( (uint64_t)(dimerTmp & 3) << (BLOCK-2)); 
			dimerTmp >>= 2; ++dtmpCount;
			if(dtmpCount==32) dimerTmp = rcwindl;
		}
		free(dimerHash);//delete [] dimerHash;
		// end checking for DIMER
		
		// start checking for HAIRPIN
		
		uint64_t hpWin=0, hpTmp=0, hpSeed=0, hpKey=0;
		hpWin = dcwindr;
		hpTmp = dcwindm;
		hpSeed = ((uint64_t)1 << 2*hpStem) - 1;
		uint64_t *hpHash = (uint64_t *) malloc ((OLIG_LENGTH-hpStem+1)*sizeof(uint64_t));
		dtmpCount=0;
		for (int i=0; i <= OLIG_LENGTH-hpStem; ++i){
			hpKey = hpWin & hpSeed;
			hpHash[i] = hpKey;
			hpWin = (hpWin >> 2) | ( (uint64_t)(hpTmp & 3) << (BLOCK-2)); 
			hpTmp >>= 2; ++dtmpCount;
			if(dtmpCount==32) hpTmp = dcwindl;
		}
		//for (int i=0; i <= OLIG_LENGTH-hpStem; ++i)
		//printf("hpHash[%d]: %16llX \n",i, hpHash[i]);
		hpWin = rcwindr;
		hpTmp = rcwindm;
		for (int i=0; i <= OLIG_LENGTH-hpStem; ++i){
			hpKey = hpWin & hpSeed;
			for (int j=0; j <= OLIG_LENGTH-hpStem; ++j)
				if (hpKey==hpHash[j] && abs(i-j)>=hpminLoop && abs(i-j)<=hpmaxLoop){
					//printf("Hairpin found!\n");
					//printf("dcwind: %16llX \t %16llX \t %16llX %16llX \n", dcwindl, dcwindm, dcwindr, hpHash[j]);
					//printf("rcwind: %16llX \t %16llX \t %16llX %16llX \n", rcwindl, rcwindm, rcwindr, hpKey);
					free(hpHash);//delete [] hpHash;
					return false;
				}
			hpWin = (hpWin >> 2) | ( (uint64_t)(hpTmp & 3) << (BLOCK-2)); 
			hpTmp >>= 2; ++dtmpCount;
			if(dtmpCount==32) hpTmp = rcwindl;
		}
		free(hpHash);//delete [] hpHash;
		// end checking for HAIRPIN
		
	}
	else{
		ehitr = get64bit(ebyte, gCoded); ehitl = get64bit(ebyte - 8, gCoded);
		ehitr >>= 2*(3-eoff);
		ehitr |= (ehitl & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		ehitl >>= 2*(3-eoff);// Remember ehitl is not a correct 32Mer (left side). it needs eoff char from ehitl2
		if (OLIG_LENGTH > 64 - (3-eoff) && ebyte-16 >=0){
			uint8_t elbyte = gCoded[ebyte-16];
			ehitl |= ((uint64_t)elbyte & (uint64_t) ( (1 << 2*(3-eoff)) - 1)) << (BLOCK - 2*(3-eoff));
		}
		
		dcwindr = ehitr;
		tpwind = ~dcwindr;
		for(int i=1; i <= 32; ++i){
			rcwindl = (rcwindl << 2) | (tpwind & 3);
			tpwind >>= 2;	
		}
		
		tpwind = 0;
		dcwindl = ehitl;
		tpwind = ~dcwindl;
		for(int i=1; i <= 18; ++i){
			rcwindr = (rcwindr << 2) | (tpwind & 3);
			tpwind >>= 2;	
		}		
		rcwindr |= ((uint64_t)(rcwindl & (((uint64_t)1 << 2*(32-(OLIG_LENGTH % 32)))-1)) << 2*(OLIG_LENGTH % 32));
		rcwindl >>= 2*(32-(OLIG_LENGTH % 32));
		
		dcwindl &= (( (uint64_t)1 << 2*(OLIG_LENGTH % 32))-1);	
		//printf("dcwind: %16llX \t %16llX \n", dcwindl, dcwindr);
		//printf("rcwind: %16llX \t %16llX \n", rcwindl, rcwindr);
		
		// start checking for DIMER
		uint64_t dimerWin=0, dimerTmp=0, dimerSeed=0, dimerKey=0;
		dimerWin = dcwindr;
		dimerTmp = dcwindl;
		dimerSeed = ((uint64_t)1 << 2*dimerLen) - 1;
		//printf("dimerWin: %16llX \n", dimerWin);
		uint64_t *dimerHash = (uint64_t *) malloc((OLIG_LENGTH-dimerLen+1)*sizeof(uint64_t));
		for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i){
			dimerKey = dimerWin & dimerSeed;
			dimerHash[i] = dimerKey;
			dimerWin = (dimerWin >> 2) | ( (uint64_t)(dimerTmp & 3) << (BLOCK-2)); 
			dimerTmp >>= 2;
		}
		/*for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i)
		 printf("dimerHash[%d]: %16llX \n",i, dimerHash[i]);*/
		dimerWin = rcwindr;
		dimerTmp = rcwindl;
		for (int i=0; i <= OLIG_LENGTH-dimerLen; ++i){
			dimerKey = dimerWin & dimerSeed;
			//printf("dimerKey: %16llX \n", dimerKey);
			// start core
			for (int j=0; j <= OLIG_LENGTH-dimerLen; ++j){
				int strGncy = 0;
				uint64_t xnorWins = ~(dimerKey ^ dimerHash[j]);
				xnorWins &= dimerSeed;
				for (int k=0; k < (dimerLen-1)/8+1; ++k)
					strGncy += preNOR[((uint16_t) (xnorWins >> 16*k))];
				//cout << (float)strGncy / (float)(dimerLen) << "\n";
				if ( (float)strGncy / (float)(dimerLen) > (dimerStg*1.0)/100.0) {
					//printf("secondary structure found!\n");
					//printf("dcwind: %16llX \t %16llX \t %16llX \n", dcwindl, dcwindr, dimerHash[j]);
					//printf("rcwind: %16llX \t %16llX \t %16llX \n", rcwindl, rcwindr, dimerKey);
					free(dimerHash);//delete [] dimerHash;
					return false;
				}
			}
			// end core
			dimerWin = (dimerWin >> 2) | ( (uint64_t)(dimerTmp & 3) << (BLOCK-2)); 
			dimerTmp >>= 2;
		}
		free(dimerHash);//delete [] dimerHash;
		// end checking for DIMER
		
		// start checking for HAIRPIN
		uint64_t hpWin=0, hpTmp=0, hpSeed=0, hpKey=0;
		hpWin = dcwindr;
		hpTmp = dcwindl;
		hpSeed = ((uint64_t)1 << 2*hpStem) - 1;
		//printf("dimerWin: %16llX \n", dimerWin);
		uint64_t *hpHash = (uint64_t *) malloc ((OLIG_LENGTH-hpStem+1)*sizeof(uint64_t));
		for (int i=0; i <= OLIG_LENGTH-hpStem; ++i){
			hpKey = hpWin & hpSeed;
			hpHash[i] = hpKey;
			hpWin = (hpWin >> 2) | ( (uint64_t)(hpTmp & 3) << (BLOCK-2)); 
			hpTmp >>= 2;
		}
		hpWin = rcwindr;
		hpTmp = rcwindl;
		for (int i=0; i <= OLIG_LENGTH-hpStem; ++i){
			hpKey = hpWin & hpSeed;
			for (int j=0; j <= OLIG_LENGTH-hpStem; ++j)
				if (hpKey==hpHash[j] && abs(i-j)>=hpminLoop && abs(i-j)<=hpmaxLoop){
					free(hpHash);//delete [] hpHash;
					return false;
				}
			hpWin = (hpWin >> 2) | ( (uint64_t)(hpTmp & 3) << (BLOCK-2)); 
			hpTmp >>= 2;
		}
		free(hpHash);//delete [] hpHash;
		// end checking for HAIRPIN		
		
	}
	return true;
}

long CalcPolyNCorrection(int gene, long pos, int NcorSize, long* NcorPos, int* Ncor, int* geneNcor){
    if(debugLevel > 3) {fprintf(stderr,"Call to CalcPolyNCorrection at %ld\n",pos);}
    int RefIndex = -1;
    int l = 0;
    int r = NcorSize - 1;
    while(l <= r){
        int m = (l + r) / 2;
        if(NcorPos[m] > pos){
            r = m - 1;
        } else if(NcorPos[m] < pos){
            l = m + 1;
            RefIndex = m;
        } else if(NcorPos[m] == pos){
            l = r+1;
            RefIndex = m;
            //This shouldn't happen: pos is the last bp in an oligo. Every
            //element of NcorPos is the first standard base after a poly-N, no
            //valid oligos contain a non-standard base.
            //Therefore if pos is in NcorPos, then the oligo defined by pos is
            //invalid, and shouldn't be printed, and shouldn't be calling this
            //function.
            printf("WARNING: Invalid oligo in sequence %d\n",gene);
        }
    }
    return (RefIndex >= 0) ? (Ncor[RefIndex] - geneNcor[gene]) : 0;
}

long DistanceToNearestN(int gene, long pos, int NcorSize, long* NcorPos, int* Ncor, int* geneNcor, bool bUpstream){
    if(debugLevel > 3) {fprintf(stderr,"Call to DistanceToNearestN %ld\n",pos);}
    int RefIndex = -1;
    int l = 0;
    int r = NcorSize - 1;
    while(l <= r){
        int m = (l + r) / 2;
        if(NcorPos[m] > pos){
            r = m - 1;
            if(!bUpstream){
                RefIndex = m;
            }
        } else if(NcorPos[m] < pos){
            l = m + 1;
            if(bUpstream){
                RefIndex = m;
            }
        } else if(NcorPos[m] == pos){
            l = r+1;
            RefIndex = m;
            //This shouldn't happen: pos is the last bp in an oligo. Every
            //element of NcorPos is the first standard base after a poly-N, no
            //valid oligos contain a non-standard base.
            //Therefore if pos is in NcorPos, then the oligo defined by pos is
            //invalid, and shouldn't be printed, and shouldn't be calling this
            //function.
            printf("WARNING: Invalid oligo in sequence %d\n",gene);
        }
    }
    return (RefIndex >= 0) ? abs(pos - NcorPos[RefIndex]) - 1 : gSize + 1; 
}

void print_Oligo(OLIGO olig, FILE* output){
    if(debugLevel > 2) {fprintf(stderr,"Call to print_Oligo(struct) at %ld in sequence %d\n",olig.uDist,olig.gene);}
    unsigned char* uPadSeq = olig.sequence;
    unsigned char* oligSeq = olig.sequence + olig.uPad;
    unsigned char* dPadSeq = oligSeq + olig.length;
    fprintf(output, "Target sequence:%d\t%.*s\tLength:%d\tTm:%0.2f\tDistance(from 5'):%d\tDistance(from 3'):%d",olig.gene+1,olig.length,oligSeq,olig.length,olig.Tm,olig.uDist,olig.dDist);
    if(olig.uPad > 0) fprintf(output,"\tPadding(5' %dbp):%.*s",olig.uPad,olig.uPad,uPadSeq);
    if(olig.dPad > 0) fprintf(output,"\tPadding(3' %dbp):%.*s",olig.dPad,olig.dPad,dPadSeq);
    fprintf(output,"\n");
}

void initOLIGO(OLIGO* self, int gene, long pos, int nOligTiles, long* geneLoc, unsigned char* gSequence, float* tmarray, long curNcor, long uNDist, long dNDist){
    int tiledLength = tilingLen + nOligTiles-1;
    self->gene = gene;
    self->length = tiledLength;
    float avgTm = 0;
    for(long i = pos - nOligTiles + 1; i <= pos; ++i){
        avgTm += tmarray[i-gseqOff] * 1.0005;
    }
    self->Tm = avgTm / nOligTiles / 1.0005;
    self->uDist = pos - gseqOff-geneLoc[gene] - tiledLength + 1 + curNcor;
    self->dDist = geneLoc[gene+1]-(pos-gseqOff)- 1 + curNcor;
    if(gene == 0) {
        self->uDist += gseqOff;
        self->dDist += gseqOff;
    }
    if(paddingLen > 0){
        self->uPad = (paddingLen < self->uDist) ? paddingLen : self->uDist;
        self->uPad = (self->uPad < uNDist) ? self->uPad : uNDist;
        self->dPad = (paddingLen < self->dDist) ? paddingLen : self->dDist;
        self->dPad = (self->dPad < dNDist) ? self->dPad : dNDist;
    } else {
        self->uPad = 0;
        self->dPad = 0;
    }
    self->sequence = gSequence + (pos-tiledLength + 1 - self->uPad) * sizeof(unsigned char);
    //if(gene == 22 && self->dDist == 3) printf("dPad: %d\n",self->dPad);
}

void print_Oligo(int gene, long pos, int nOligTiles, long* geneLoc, int* geneLen, FILE* output, unsigned char* gSequence, float* tmarray,int NcorSize, long* NcorPos, int* Ncor, int* geneNcor)
{
    if(debugLevel > 2) {fprintf(stderr,"Call to print_Oligo at %ld\n",pos);}
    long curNcor = CalcPolyNCorrection(gene,pos,NcorSize,NcorPos,Ncor,geneNcor);
    int tiledLength = tilingLen + nOligTiles-1;
    //if(tiledLength > OLIG_LENGTH){
    //    printf("Gene: %d, pos: %d, tiles: %d\n",gene,pos,nOligTiles);
    //}
    fprintf(output, "_Target sequence:%d\t",gene+1);
    for (long i = pos - tiledLength + 1; i <= pos; ++i){
    	fprintf(output, "%c", gSequence[i]);
    }
    fprintf(output, "\tLength:%d", tiledLength); 
    float avgTm = 0;
    for(long i = pos - nOligTiles + 1; i <= pos; i++){
        avgTm += tmarray[i-gseqOff];
    }
    avgTm /= nOligTiles;
    fprintf(output, "\tTm%s:%.2f",(nOligTiles > 1) ? "(avg)" : "", avgTm);
    fprintf(output, "\tDistance(from 5'):%ld",
        pos-gseqOff-geneLoc[gene] - tiledLength + 1 + curNcor);
    fprintf(output, "\tDistance(from 3'):%ld",
        geneLoc[gene+1]-(pos-gseqOff)-1 + curNcor);
    if(paddingLen > 0){
        long LPad = pos - tiledLength - paddingLen + 1;
        if(LPad < geneLoc[gene]) LPad = geneLoc[gene];
        fprintf(output,"\tPadding(5' %dbp):",pos - LPad - tiledLength+1);
        for (long i = LPad; i <= pos - tiledLength; ++i){
            fprintf(output, "%c", gSequence[i]);
        }
        long RPad = pos + paddingLen;
        if(RPad >= geneLoc[gene + 1]) RPad = geneLoc[gene+1] - 1;
        fprintf(output,"\tPadding(3' %dbp):",RPad - pos);
        for (long i = pos+1; i <= RPad; ++i){
            fprintf(output, "%c", gSequence[i]);
        }
    }
    fprintf(output, "\n");
    //++olig_total;
}

void print_Intervals(long* geneLoc, int* geneLen, FILE* output, unsigned char* gSequence, float* tmarray, uint8_t* simarray,int NcorSize,long* NcorPos,int* Ncor,int* geneNcor){
    if(debugLevel > 1) {fprintf(stderr,"Call to print_Intervals\n");}
    int curSeq = 0;
    long pos = OLIG_LENGTH -1;
    while(pos < gSize){
        if(simarray[pos] == 0){
            long curNcor = CalcPolyNCorrection(curSeq,pos,NcorSize,NcorPos,Ncor,geneNcor);
            fprintf(output, "Target sequence:%d \t",curSeq+1);
            long i = pos - OLIG_LENGTH + 1;
            while((i <= pos || simarray[i] == 0) && i < gSize){
                fprintf(output, "%c", gSequence[i]);
                i++;
            }
            int regLen = OLIG_LENGTH + (i-pos-1);
            fprintf(output, "\tLength:%d", regLen); 
            fprintf(output, "\tTm:%.2f(first)", tmarray[pos-gseqOff]); 
            fprintf(output, "\tDistance(from 5'):%ld",
                pos-gseqOff-geneLoc[curSeq]-OLIG_LENGTH + 1 + curNcor);
            fprintf(output, "\tDistance(from 3'):%ld",
            geneLoc[curSeq+1]-(i-gseqOff-1)-1,+curNcor);
            long LPad = pos - OLIG_LENGTH - paddingLen + 1;
            if(LPad < geneLoc[curSeq]) LPad = geneLoc[curSeq];
            int LPadLen = pos - LPad - OLIG_LENGTH +1;
            if(LPadLen){
                fprintf(output,"\tPadding(5' %dbp):",LPadLen);
                for (long j = LPad; j <= pos - OLIG_LENGTH; ++j){
                    fprintf(output, "%c", gSequence[j]);
                }
            }
            if(i > pos) pos = i-1;
            long RPad = pos + paddingLen;
            if(RPad >= geneLoc[curSeq + 1]) RPad = geneLoc[curSeq+1] - 1;
            int RPadLen = RPad - pos;
            if(RPadLen){
                fprintf(output,"\tPadding(3' %dbp):",RPadLen);
                for (long j = pos+1; j <= RPad; ++j){
                    fprintf(output, "%c", gSequence[j]);
                }
            }
            fprintf(output, "\n");
            ++olig_total;
        }
        if(pos >= geneLoc[curSeq+1]){++curSeq;}
        pos++;
    }
}
long maxSize = 0;
unsigned long intensive_Homology(uint8_t *gCoded, uint8_t *simarray, unsigned char *gSequence, long *geneLoc, float *tmarray, int *geneLen,int NcorSize, long* NcorPos, int* Ncor, int* geneNcor, OLIGO* &OLIGList){
        if(debugLevel > 0) {fprintf(stderr,"Call to intensive_Homology\n");}
	unsigned long memoryUsed=0;
        time_t timer_s, timer_f;
	int nthreads, tid;
	if (MAX_NUM_THREAD >= 8)	
		omp_set_num_threads(8);
	else 
		omp_set_num_threads(MAX_NUM_THREAD);

        entry* hashTableList[8];
        uint64_t seedSeqList[8];
        uint64_t seedLenList[8];
        int* oligsperGene = NULL;
        int* geneSize = NULL;
        OLIGO** Gene_OligList = NULL;
        oligsperGene = (int*) malloc(geneno * sizeof(int));
        Gene_OligList = (OLIGO**) malloc(geneno * sizeof(OLIGO*));
        memoryUsed += geneno * sizeof(OLIGO*);
        for(int i = 0; i < geneno; ++i){
            oligsperGene[i] = 0;
        }
        if(maxOligs == -1){
            geneSize = (int*) malloc(geneno * sizeof(int));
            for(int i = 0; i < geneno; i++){
                geneSize[i] = 1;
                Gene_OligList[i] = (OLIGO*) malloc(sizeof(OLIGO));
            }
            memoryUsed += geneno * (sizeof(OLIGO) + sizeof(int));

        } else {
            for(int i = 0; i < geneno; i++){
                Gene_OligList[i] = (OLIGO*) malloc(maxOligs * sizeof(OLIGO));
            }
            memoryUsed += maxOligs * geneno * sizeof(OLIGO);
        }
        timer_s = time(NULL);
        int chunk = 1;
        int curTable = 0;
#pragma omp parallel for shared(gCoded,hashTableList,seedSeqList,seedLenList) private(curTable) schedule(dynamic,chunk)
        for(curTable = 0; curTable < 8; ++curTable){
            hashTableList[curTable] = (entry*) malloc(PRIME_IH * sizeof(entry));
            seedSeqList[curTable] = get64seed(seed_Phase2[curTable]);
            seedLenList[curTable] = strlen(seed_Phase2[curTable]);
            create_Hash(hashTableList[curTable],gCoded, seedSeqList[curTable],seedLenList[curTable],nbytes,PRIME_IH);
        }

        memoryUsed += 8*PRIME_IH * sizeof(entry);

        unsigned long sumPos = 0;
	for (long i=0; i < PRIME_IH; ++i) {
                for(int j = 0; j < 8; ++j)
		sumPos += hashTableList[j][i].size;
	
        }
	memoryUsed = sumPos * sizeof (long) + 8;

        timer_f = time(NULL);
	printf("Intensive Homology Hash Tables Created in %ds\n",(timer_f - timer_s));
        //Hash Table Structure
        //An array of size 8, one for each Key, each is composed of
        //An array structs of size PRIME_IH
        //A Key, and size, which desribes the position array
        //for(int i=0; i<8;i++){
        //    print_HashTable(hashTableList[i]);
        //}
        
        uint8_t** tmp_simarray = (uint8_t**) malloc(sizeof(uint8_t*) * MAX_NUM_THREAD);
        for(int i = 0; i < MAX_NUM_THREAD; i++){
            tmp_simarray[i] = (uint8_t*) malloc(sizeof(uint8_t) * gSize);
        }
        memoryUsed += MAX_NUM_THREAD * gSize * sizeof(uint8_t);

	olig_total = 0;
	chunk = (geneno < MAX_NUM_THREAD) ? 1 : geneno / MAX_NUM_THREAD;

        //printf("PreSearch\n");

	if (MAX_NUM_THREAD >= geneno)	
	    omp_set_num_threads(geneno);
	else 
	    omp_set_num_threads(MAX_NUM_THREAD);
        int i;
#pragma omp parallel for shared(simarray) private(i) schedule(dynamic,chunk)
	for (i=0; i < geneno; ++i){
            //printf("gene: %d\t",i);
            #pragma omp flush(simarray)
            {
                memcpy(tmp_simarray[i % MAX_NUM_THREAD],simarray,gSize);
            }
            int seenZero;
            long bestZero, seenPlace, bestPlace, canPlace;
            long match_pos;
            long lastPlace = -1;
            int count = 0;

	    do{
                match_pos = -1;
		seenZero = 0, bestZero = 0, seenPlace = -1, bestPlace = -1, canPlace=-1;
                if(maxOligs == -1){ //Find the first potential oligo
                    long j = geneLoc[i] + OLIG_LENGTH -1;
                    while(tmp_simarray[i %MAX_NUM_THREAD][j] != 0 & j < geneLoc[i+1]) j++;
                    if(j < geneLoc[i+1]){
                        bestPlace = j;
                    }
                } else { // Find the oligo in the middle of the longest still valid stretch
		    for(long j = geneLoc[i] + OLIG_LENGTH - 1; j < geneLoc[i+1]; ++j){
		        if (tmp_simarray[i % MAX_NUM_THREAD][j] == 0){
		        	++seenZero;
		        	seenPlace = j;
		        	if (seenZero > bestZero) {bestZero = seenZero; bestPlace = seenPlace;}
		        }
		        else 
		        	seenZero = 0;
		    }
                }
		if (bestPlace != -1){
		    canPlace = (tilingLen == OLIG_LENGTH) ? bestPlace - bestZero/2 : bestPlace;
                    //if(lastPlace == canPlace){
                    //    printf("cp: %ld\n",canPlace);
                    //}
                    lastPlace = canPlace;
		    if (secStr){
                        bool res = secondary_Structure(gCoded, canPlace);
                        if(!res) match_pos = -3; // Invalid secondary_Structure result 
                    }
                    for (int j = 0; j < 8  && match_pos == -1; ++j)
                    {
                        bool res = second_Check(seedSeqList[j], seedLenList[j], hashTableList[j], gCoded, canPlace,geneLoc[i],geneLoc[i+1]-1, tmp_simarray[i % MAX_NUM_THREAD]);
                        if(!res) match_pos = 0;
                    }
                    if (match_pos == -1){
                        long max_r, cur_r;
                        int max_run = -1;
                        int cur_run = 0;
                        for(long j = 0; j < nTiles; ++j){
                            
                            if(!tmp_simarray[i % MAX_NUM_THREAD][canPlace-j]){
                                if(cur_run == 0){
                                    cur_r = canPlace -j;
                                }
                                cur_run++;
                                if(cur_run > max_run){
                                    max_run = cur_run;
                                    max_r = cur_r;
                                }
                            } else {
                                cur_run = 0;
                            }
                        }
                        canPlace = max_r;
                        for(int j = 0; j < max_run; ++j){
                            tmp_simarray[i % MAX_NUM_THREAD][canPlace-j] = 1;
                        }
                        if(maxOligs == -1 && oligsperGene[i] >= geneSize[i]){
                            geneSize[i] *= 2;
                            Gene_OligList[i] = (OLIGO*) realloc(Gene_OligList[i],geneSize[i] * sizeof(OLIGO));
                        }
                        long gPos = canPlace+gseqOff;
                        long Ncorval = CalcPolyNCorrection(i,gPos,NcorSize,NcorPos,Ncor,geneNcor);
                        long uNDist = DistanceToNearestN(i,gPos,NcorSize,NcorPos,Ncor,geneNcor,1);
                        long dNDist = DistanceToNearestN(i,gPos,NcorSize,NcorPos,Ncor,geneNcor,0);
                        initOLIGO(&Gene_OligList[i][oligsperGene[i]++],i,gPos,max_run,geneLoc,gSequence,tmarray,Ncorval,uNDist,dNDist);
		    }else if(match_pos == -3){
                        tmp_simarray[i % MAX_NUM_THREAD][canPlace] = 1;
                        match_pos = -1;
		    }
		}else 
		    match_pos = -2; //All positions invalid

	    }while(match_pos != -2 && (maxOligs == -1 || oligsperGene[i] < maxOligs)); //While matches are found/possible
            #pragma omp flush(simarray)
            for(long j = 0; j < gSize; j++){
                simarray[j] |= tmp_simarray[i % MAX_NUM_THREAD][j];
            }
	}
        //printf("PostSearch\n");
        for(int i = 0; i < MAX_NUM_THREAD; i++){
            free(tmp_simarray[i]);
        }
        free(tmp_simarray);
        //printf("PostFree\n");
        olig_total = 0;
        for(int i = 0; i < geneno; i++){
            olig_total += oligsperGene[i];
        }
        //printf("PostCount\n");
        OLIGList = (OLIGO*) malloc(olig_total * sizeof(OLIGO));
        memoryUsed += olig_total * sizeof(OLIGO);
        int curOlig = 0;
        //printf("PreLoop\n");
        for(int i = 0; i < geneno; i++){
            memoryUsed += oligsperGene[i] * sizeof(OLIGO);
            for(int j = 0; j < oligsperGene[i]; j++){
                //printf("Olig: %d, Gene: %d, OligInGene: %d\n",curOlig,i,j);
                memcpy(&OLIGList[curOlig++],&Gene_OligList[i][j],sizeof(OLIGO));
            }
            free(Gene_OligList[i]);
        }
        free(Gene_OligList);
        //printf("PostCopy\n");
	for(int i = 0; i < 8; ++i){
            for(int j = 0; j < PRIME_IH; j++){
                free(hashTableList[i][j].pos);
            }
            free(hashTableList[i]);
        }
        free(oligsperGene);
        free(geneSize);
		
	return memoryUsed;
}

void pre_Process(){
        if(debugLevel > 1) {fprintf(stderr,"Call to pre_Process\n");}
	for (int i=0; i < 65536;++i)
		preNOR[i]=0;
	for (int i=0; i < 65536;++i){
		uint16_t tmp=i;
		for (int pass=0; pass<8;++pass){
			if ((tmp & 3) == 3) 
				++preNOR[i];
			tmp >>= 2;
		}
	}
}

unsigned long pickPrime(){
        if(debugLevel > 1) {fprintf(stderr,"Call to pickPrime\n");}
	if(gSize <  30000000)
		return (34040383);
	else if (gSize < 113493747)
		return (113493747);
	const unsigned long primeTable[] = {
		1114523, 1180043, 1245227, 1310759, 1376447, 1442087, 1507379,
		1573667, 1638899, 1704023, 1769627, 1835027, 1900667, 1966127,
		2031839, 2228483, 2359559, 2490707, 2621447, 2752679, 2883767,
		3015527, 3145739, 3277283, 3408323, 3539267, 3670259, 3801143,
		3932483, 4063559, 4456643, 4718699, 4980827, 5243003, 5505239,
		5767187, 6029603, 6291563, 6553979, 6816527, 7079159, 7340639,
		7602359, 7864799, 8126747, 8913119, 9437399, 9962207, 10485767,
		11010383, 11534819, 12059123, 12583007, 13107923, 13631819, 14156543,
		14680067, 15204467, 15729647, 16253423, 17825999, 18874379, 19923227,
		20971799, 22020227, 23069447, 24117683, 25166423, 26214743, 27264047,
		28312007, 29360147, 30410483, 31457627, 32505983, 35651783, 37749983,
		39845987, 41943347, 44040383, 46137887, 48234623, 50331707, 52429067,
		54526019, 56623367, 58720307, 60817763, 62915459, 65012279, 71303567,
		75497999, 79691867, 83886983, 88080527, 92275307, 96470447, 100663439,
		104858387, 109052183, 113246699, 117440699, 121635467, 125829239,
		130023683, 142606379, 150994979, 159383759, 167772239, 176160779,
		184549559, 192938003, 201327359, 209715719, 218104427, 226493747,
		234882239, 243269639, 251659139, 260047367, 285215507, 301989959,
		318767927, 335544323, 352321643, 369100463, 385876703, 402654059,
		419432243, 436208447, 452986103, 469762067, 486539519, 503316623,
		520094747, 570425399, 603979919, 637534763, 671089283, 704643287,
		738198347, 771752363, 805307963, 838861103, 872415239, 905971007,
		939525143, 973079279, 1006633283, 1040187419, 1140852767, 1207960679,
		1275069143, 1342177379, 1409288183, 1476395699, 1543504343, 1610613119,
		1677721667, 1744830587, 1811940419, 1879049087, 1946157419, 2013265967,
		2080375127, 2281701827, 2415920939, 2550137039, 2684355383, 2818572539,
		2952791147, 3087008663, 3221226167, 3355444187, 3489661079, 3623878823,
		3758096939, 3892314659, 4026532187, 4160749883, 4563403379, 4831838783,
		5100273923, 5368709219, 5637144743, 5905580687, 6174015503, 6442452119,
		6710886467, 6979322123, 7247758307, 7516193123, 7784629079, 8053065599,
		8321499203, 9126806147, 9663676523, 10200548819, 10737418883,
		11274289319, 11811160139, 12348031523, 12884902223, 13421772839,
		13958645543, 14495515943, 15032386163, 15569257247, 16106127887,
		16642998803, 18253612127, 19327353083, 20401094843, 21474837719,
		22548578579, 23622320927, 24696062387, 25769803799, 26843546243,
		27917287907, 28991030759, 30064772327, 31138513067, 32212254947,
		33285996803, 36507222923, 38654706323, 40802189423, 42949673423,
		45097157927, 47244640319, 49392124247, 51539607599, 53687092307,
		55834576979, 57982058579, 60129542339, 62277026327, 64424509847,
		66571993199, 73014444299, 77309412407, 81604379243, 85899346727,
		90194314103, 94489281203, 98784255863, 103079215439, 107374183703,
		111669150239, 115964117999, 120259085183, 124554051983, 128849019059,
		133143986399, 146028888179, 154618823603, 163208757527, 171798693719,
		180388628579, 188978561207, 197568495647, 206158430447, 214748365067,
		223338303719, 231928234787, 240518168603, 249108103547, 257698038539,
		266287975727, 292057776239, 309237645803, 326417515547, 343597385507,
		360777253763, 377957124803, 395136991499, 412316861267, 429496730879,
		446676599987, 463856468987, 481036337207, 498216206387, 515396078039,
		532575944723, 584115552323, 618475290887, 652835029643, 687194768879,
		721554506879, 755914244627, 790273985219, 824633721383, 858993459587,
		893353198763, 927712936643, 962072674643, 996432414899, 1030792152539,
		1065151889507, 1168231105859, 1236950582039, 1305670059983, 1374389535587,
		1443109012607, 1511828491883, 1580547965639, 1649267441747, 1717986918839,
		1786706397767, 1855425872459, 1924145348627, 1992864827099, 2061584304323,
		2130303780503, 2336462210183, 2473901164367, 2611340118887, 2748779070239,
		2886218024939, 3023656976507, 3161095931639, 3298534883999, 3435973836983,
		3573412791647, 3710851743923, 3848290698467, 3985729653707, 4123168604483,
		4260607557707, 4672924419707, 4947802331663, 5222680234139, 5497558138979,
		5772436047947, 6047313952943, 6322191860339, 6597069767699, 6871947674003,
		7146825580703, 7421703488567, 7696581395627, 7971459304163, 8246337210659,
		8521215117407, 9345848837267, 9895604651243, 10445360463947,
		10995116279639, 11544872100683, 12094627906847, 12644383722779,
		13194139536659, 13743895350023, 14293651161443, 14843406975659,
		15393162789503, 15942918604343, 16492674420863, 17042430234443,
		18691697672867, 19791209300867, 20890720927823, 21990232555703,
		23089744183799, 24189255814847, 25288767440099, 26388279068903,
		27487790694887, 28587302323787, 29686813951463, 30786325577867,
		31885837205567, 32985348833687, 34084860462083, 37383395344739,
		39582418600883, 41781441856823, 43980465111383, 46179488367203,
		48378511622303, 50577534878987, 52776558134423, 54975581392583,
		57174604644503, 59373627900407, 61572651156383, 63771674412287,
		65970697666967, 68169720924167, 74766790688867, 79164837200927,
		83562883712027, 87960930223163, 92358976733483, 96757023247427,
		101155069756823, 105553116266999, 109951162779203, 114349209290003,
		118747255800179, 123145302311783, 127543348823027, 131941395333479,
		136339441846019, 149533581378263, 158329674402959, 167125767424739,
		175921860444599, 184717953466703, 193514046490343, 202310139514283,
		211106232536699, 219902325558107, 228698418578879, 237494511600287,
		246290604623279, 255086697645023, 263882790666959, 272678883689987,
		299067162755363, 316659348799919, 334251534845303, 351843720890723,
		369435906934019, 387028092977819, 404620279022447, 422212465067447,
		439804651111103, 457396837157483, 474989023199423, 492581209246163,
		510173395291199, 527765581341227, 545357767379483, 598134325510343,
		633318697599023, 668503069688723, 703687441776707, 738871813866287,
		774056185954967, 809240558043419, 844424930134187, 879609302222207,
		914793674313899, 949978046398607, 985162418489267, 1020346790579903,
		1055531162666507, 1090715534754863
	};
	int i=0;
	while (primeTable[i] < 2*gSize) i++;//printf("Current Prime Index: %d\n",i++);
	return primeTable[i];
}

//UNUSED
//Assumes complement is already allocated enough memory to containin
//the reverse complement of string
void append_RevComp(int stringlength, unsigned char* string, unsigned char* complement){
    if(debugLevel > 3) {fprintf(stderr,"Call to append_RevComp\n");}
    for(int i = stringlength - 1; i >= 0; --i){
       switch (string[i]){
            case 'A':
               complement[stringlength - i - 1] = 'T';
               break;
            case 'C':
               complement[stringlength - i - 1] = 'G';
               break;
            case 'G':
               complement[stringlength - i - 1] = 'C';
               break;
            case 'T':
               complement[stringlength - i - 1] = 'A';
               break;
            default:
               printf("ERROR in revcomp at position %d\n\n",i);
                exit(EXIT_FAILURE);
               break;
       }
    }
}

uint8_t* encodeSequence(uint8_t* &Coded, long sSize, unsigned char* Sequence){
    long nbytes = sSize/4;
    if(!Coded) return 0;
    uint8_t temp = 0;
    int stpoint = sSize % 4;
    for(long i = stpoint; i < sSize; i += 4) {
        for(int j = 0; j < 4; j++) 
            switch (Sequence[j+i]) {
                case 'A':
                	temp = (uint8_t) (temp << 2 | 0);
                	break;
                case 'C':
                	temp = (uint8_t) (temp << 2 | 1);
                	break;
                case 'G':
                	temp = (uint8_t) (temp << 2 | 2);
                	break;
                case 'T':
                	temp = (uint8_t) (temp << 2 | 3);
                	break;
                default:
                	cout << "Error! " << Sequence[j+i] << "is not a valid nucleotide.\n";
                	temp = (uint8_t) (temp << 2 | 0);					
            }	
        Coded[(i-stpoint)/4] = temp;
    }
    return Coded;
}



unsigned long ProcessBlacklistSeq(long sSize, unsigned char* Sequence, OLIGO* &OLIGList){
    if(olig_total < 1) return 0; // No need to do anything if there are no oligos

    unsigned long memoryEval;
    unsigned long memoryUsed = 0;
    unsigned long peakMemory = 0;

    ///// Clean Up the Sequence
    unsigned char* bSequence = (unsigned char*) malloc(sSize * sizeof(unsigned char));
    peakMemory += sizeof(unsigned char) * sSize;
    long i = 0;
    long bSize = 0;
    bool prevNonStandard = false;
    while (i < sSize){
        switch (Sequence[i]) {
    	    case 'a': case 'A':
        	bSequence[bSize++] = 'A';
                prevNonStandard = false;
        	break;
            case 'c': case 'C':
        	bSequence[bSize++] = 'C';
                prevNonStandard = false;
        	break;
            case 'g': case 'G':
        	bSequence[bSize++] = 'G';
                prevNonStandard = false;
        	break;
            case 't': case 'T':
        	bSequence[bSize++] = 'T';
                prevNonStandard = false;
        	break;
            case '\n': case '\r':
        	break;
            default: //Non-Standard nucleotides
                if(!prevNonStandard){
	    	    bSequence[bSize++] = 'A'; 
	    	    prevNonStandard = true;
                }
        }
        ++i;
    }
    bSequence = (unsigned char*) realloc(bSequence,bSize * sizeof(unsigned char));
    memoryUsed = sizeof(unsigned char) * bSize;

    ///// Encode The Sequence
    uint8_t* bCoded = (uint8_t *) malloc (bSize/4 * sizeof(uint8_t));
    if(!encodeSequence(bCoded,bSize,bSequence)) {
        cout << "Could not encode blacklist sequence\n";
        exit(EXIT_FAILURE);
    }
    memoryUsed += (bSize - bSize % 4) * sizeof(uint8_t);
    if(memoryUsed > peakMemory) peakMemory = memoryUsed;
    free(bSequence);
    memoryUsed -= sizeof(unsigned char) * bSize;
    bSize -= bSize % 4;

    int numSeeds = 8;

    //Create Hash Table
    int nthreads, tid;
    if (MAX_NUM_THREAD >= numSeeds)	
	omp_set_num_threads(numSeeds);
    else 
    	omp_set_num_threads(MAX_NUM_THREAD);
    
    entry* hashTableList[8];
    uint64_t seedSeqList[8];
    uint64_t seedLenList[8];
    unsigned long PRIME_BL = PRIME_IH;//5452619;
    
    int chunk = 1;
#pragma omp parallel for shared(bCoded) schedule(dynamic,chunk)
    for(int i = 0; i < numSeeds; ++i){
        hashTableList[i] = (entry*) malloc(PRIME_BL * sizeof(entry));
        seedSeqList[i] = get64seed(seed_Phase2[i]);
        seedLenList[i] = strlen(seed_Phase2[i]);
        create_Hash(hashTableList[i],bCoded, seedSeqList[i],seedLenList[i],bSize / 4,PRIME_BL);
    }
    

    memoryUsed += numSeeds*PRIME_BL * sizeof(entry);
    unsigned long sumPos = 0;
    for (long i=0; i < PRIME_BL; ++i) {
            for(int j = 0; j < numSeeds; ++j)
    	sumPos += hashTableList[j][i].size;
    
    }
    memoryUsed = sumPos * sizeof (long) + numSeeds;
    if(memoryUsed > peakMemory) peakMemory = memoryUsed;

    ///// Test Oligos
    bool* OLIGPass = (bool*) malloc(olig_total * sizeof(bool));
    for(int i = 0; i < olig_total; i++){
        OLIGPass[i] = true;
    }
    memoryUsed = olig_total * sizeof (bool);
    if(memoryUsed > peakMemory) peakMemory = memoryUsed;
    int j = 0;
    chunk = (olig_total < MAX_NUM_THREAD) ? 1 : olig_total / MAX_NUM_THREAD;
    omp_set_num_threads((olig_total < MAX_NUM_THREAD) ? olig_total : MAX_NUM_THREAD);
#pragma omp parallel for shared(bCoded,OLIGList,olig_total) private(j) schedule(dynamic,chunk) 
    for (j=0; j < olig_total; ++j){
        for (int k = 0; k < numSeeds  && OLIGPass[j]; ++k)
        {
            bool res = blacklist_Check(seedSeqList[k], seedLenList[k], hashTableList[k], bSize, bCoded, OLIGList[j]);
            if(!res) OLIGPass[j] = false;
        }
    }

    ////Move OLIGOS from the end to fill in positions of oligos which match the blacklist
    for(int i = 0; i < olig_total; i++){
        if(!OLIGPass[i]){
            if(i != olig_total-1) memcpy(&OLIGList[i],&OLIGList[olig_total-1],sizeof(OLIGO));
            OLIGPass[i] = OLIGPass[olig_total -1];
            i--;
            olig_total--;
        }
    }

    OLIGList = (OLIGO*) realloc(OLIGList,olig_total * sizeof(OLIGO));

    free(bCoded);
    free(OLIGPass);
    for(int i = 0; i < numSeeds; ++i){
        for(int j = 0; j < PRIME_BL; j++){
            free(hashTableList[i][j].pos);
        }
        free(hashTableList[i]);
    }

    return peakMemory;
}

unsigned long ProcessBlacklistFile(char * file,OLIGO* &OLIGList){
    if(olig_total < 1) return 0; //No need to do anything if there are no oligos
    FILE* bf;
    if((bf = fopen(file, "rb")) == NULL){
        printf("Error! Could not open blacklist file - %s: Skipping\n",file);
        return 0 ;
    }
    long initSize = 1000;
    unsigned long memoryEval;
    unsigned long memoryUsed = 0;
    unsigned long peakMemory = 0;
    long curSize = initSize;
    long sSize = 0;
    unsigned char* tempSequence = (unsigned char*) malloc(sizeof(unsigned char) * curSize);
    memoryUsed += sizeof(unsigned char) * curSize;
    peakMemory = memoryUsed;
    unsigned char line[100];
    while(fgets((char*)line, 100,bf)){
        if(line[0] == '>'){
            if(sSize){
                tempSequence = (unsigned char*) realloc(tempSequence,sizeof(unsigned char) * sSize);
                memoryEval = ProcessBlacklistSeq(sSize,tempSequence, OLIGList);
                if (memoryUsed  + memoryEval > peakMemory) peakMemory = memoryUsed + memoryEval;
                free(tempSequence);
                sSize = 0;
                curSize = initSize;
                tempSequence = (unsigned char*) malloc(sizeof(unsigned char) * curSize);
                memoryUsed = sizeof(unsigned char) * curSize;
            }
            line[strlen((char*)line)-1] = '\0';
            printf("Filtering against %s\n",line);
        } else {
            if(sSize + 100 >= curSize){
                curSize *= 2;
                tempSequence = (unsigned char*) realloc(tempSequence,sizeof(unsigned char) * curSize);
                memoryUsed += sizeof(unsigned char) * curSize;
            }
            strcpy((char*)&tempSequence[sSize],(char*)line);
            sSize += strlen((char*)line);
        }
    }
    if(sSize){
        memoryEval = ProcessBlacklistSeq(sSize,tempSequence, OLIGList);
        if (memoryUsed  + memoryEval > peakMemory) peakMemory = memoryUsed + memoryEval;
        free(tempSequence);
    }
    fclose(bf);
    //cout << peakMemory << "\n";
    return peakMemory;
}

int main (int argc, char * const argv[]) {
    	cout << "------------------------------BEGIN------------------------------------\n";
	print_Time();

	time_t seconds, secondf, timers, timerf;	
	seconds = time (NULL);
	
	
		
	unsigned long memoryUsed = 0, peakUsed = 0;
	
	
	pre_Process();
	
	FILE *gf, *oligos;
	unsigned char *tempSequence, *gSequence;
	uint8_t *gCoded;
	uint8_t *simarray;
	long tempSize;
        char* blackliststr = NULL;
	
	const char *gfname = argv[1];
	/* Open genome file for reading. */
	if((gf = fopen(gfname, "rb")) == NULL) {
		cout << "Error! cannot open file: "<< gfname << endl;
		exit(EXIT_FAILURE);
	}
	
	const char *oligosname = argv[2];
	if (argv[2][0]=='-'){
		cout << "Error! no output file. \n";
		exit(EXIT_FAILURE);
	}

	
	for (int i=3; i < argc; i+=2){		
		string param = argv[i];
		if (param == "-length") {OLIG_LENGTH = atoi(argv[i+1]); tilingLen = OLIG_LENGTH;}
		else if (param == "-seqSim") seqSimilarity = atoi(argv[i+1]); 
		else if (param == "-maxMatch") consecMatch = atoi(argv[i+1]); 
		else if (param == "-maxGC") GC_CONT_MAX = atoi(argv[i+1]); 
		else if (param == "-minGC") GC_CONT_MIN = atoi(argv[i+1]); 
		else if (param == "-dimerLen") dimerLen = atoi(argv[i+1]); 
		else if (param == "-dimerSim") dimerStg = atoi(argv[i+1]); 
		else if (param == "-hairpinLen") hpStem = atoi(argv[i+1]); 				
		else if (param == "-minhpLoop") hpminLoop = atoi(argv[i+1]); 
		else if (param == "-maxhpLoop") hpmaxLoop = atoi(argv[i+1]); 				
		else if (param == "-minTm") tmmin = atoi(argv[i+1]); 
		else if (param == "-maxTm") tmmax = atoi(argv[i+1]);
		else if (param == "-rangeTm") tminterval = atoi(argv[i+1]); 
		else if (param == "-oligCon") oligc = (double)(atoi(argv[i+1]))/1000000000.0; 
		else if (param == "-saltCon") saltc = (double)(atoi(argv[i+1]))/1000.0;	
		else if (param == "-secStr") secStr = true;	
                else if (param == "-maxOligs") maxOligs = atoi(argv[i+1]); 
                else if (param == "-paddingLen") paddingLen = atoi(argv[i+1]);
                else if (param == "-tilingLen") tilingLen = atoi(argv[i+1]);
                else if (param == "-maxThread") maxThread = atoi(argv[i+1]);
                else if (param == "-debug") debugLevel=atoi(argv[i+1]);
                else if (param == "-blacklist") blackliststr=argv[i+1];
		else {
			cout << "Error! in parameters!\n";
			exit(EXIT_FAILURE);
		}
	}
        if(tilingLen > OLIG_LENGTH){
            cout << "WARNING: tilingLen must be less than or equal to length, continuing with tilingLen as length\n";
            tilingLen = OLIG_LENGTH;
        }
        nTiles = (OLIG_LENGTH - tilingLen) + 1;
        //Calculate number of blocks required to store an oligo
	nBlocks = OLIG_LENGTH/(BLOCK/2);
        if(OLIG_LENGTH % (BLOCK/2) != 0) ++nBlocks;
        lchars = nBlocks * (BLOCK / 2) - OLIG_LENGTH - 1;
        nNORs = OLIG_LENGTH/8;
        NORoff = OLIG_LENGTH % 8;
        tNORs = tilingLen / 8;
        tNORoff = tilingLen % 8;
        //Calculate the number of identities that will give the
        //maximum similarity
        maxseqIdentity = seqSimilarity * OLIG_LENGTH / 100 + 1;
        maxtilingIdentity = (seqSimilarity * tilingLen / 100 + 1);
        mintilingIdentityDiff = maxtilingIdentity - maxseqIdentity;
        maxtilingIdentityDiff = OLIG_LENGTH - ((1-seqSimilarity) * tilingLen / 100 + 1) - maxseqIdentity;
		
	/* Get the genome file size. */
	if(fseek(gf, 0, SEEK_END) == 0) {
		tempSize = ftell(gf);//*2;
		rewind(gf);
		if(tempSize < 0) { 
			cout << "Error! Cannot ftell: " << gfname;
			perror(NULL);
			exit(EXIT_FAILURE);
		}
	} else {
		cout << "Error! Cannot fseek: " << gfname;
		perror(NULL);
		exit(EXIT_FAILURE);
	}
        
        //printf("Tempsize: %ld\n",tempSize);
	
	/* Allocate memory to genome sequence array. */
	tempSequence = (unsigned char *) malloc(tempSize * sizeof(unsigned char));
	memoryUsed += (tempSize) * sizeof(unsigned char);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	gSequence = (unsigned char *) malloc(tempSize * sizeof(unsigned char));
	memoryUsed += (tempSize) * sizeof(unsigned char);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	simarray = (uint8_t *) malloc (tempSize * sizeof(uint8_t));
	memoryUsed += (tempSize) * sizeof(uint8_t);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	for (long i=0; i<tempSize; ++i)
		simarray[i] = 0;
	
	if((tempSequence == NULL) || (gSequence == NULL)) {
		cout << "Error! Cannot allocate memory.";
		exit(EXIT_FAILURE);
	}
	
	/* Read tempSize bytes of data. */
	if(fread(tempSequence, sizeof(unsigned char), tempSize/*/2*/, gf) != tempSize){
		cout << "Cannot read from: "<< gfname << "\n";
		perror(NULL);
		exit(EXIT_FAILURE);
	}
	fclose(gf);

        //printf("temporary memory allocated\n");
        
	long gene_num = 0;
	int locSize = 8;
	int blanks=0;
	geneno=0;
	int errorno = 0;
	long i = 0;
	long *geneLoc = (long *) malloc(locSize * sizeof(long));
        int *geneLen = (int *) malloc(locSize * sizeof(int));
        int* geneNcor = (int*) malloc(locSize * sizeof(int));
        int NcorSize = (tempSize >= 1000) ? tempSize / 1000 : 8;
        long* NcorPos = (long*) malloc(NcorSize * sizeof(long));
        int* Ncor = (int*) malloc(NcorSize * sizeof(int));
        int curNcor = 0;
        int curNcorIndex = 0;
        bool prevNonStandard = false;
	while (i < tempSize){
	    switch (tempSequence[i]) {
	    	case '>':
                    if (geneno > 0){
                        geneLen[geneno-1] = gene_num - geneLoc[geneno - 1]; 
                    }
	    	    if (geneno >= locSize){
	    	    	locSize *= 2;
	    	    	geneLoc = (long *) realloc(geneLoc, locSize * sizeof(long));
                        geneLen = (int *) realloc(geneLen, locSize * sizeof(int));
                        geneNcor = (int *) realloc(geneNcor, locSize * sizeof(int));
	    	    }
	    	    geneLoc[geneno++] = gene_num;
                    geneNcor[geneno-1] = curNcor;
	    	    while (tempSequence[i] != '\n' && tempSequence[i] != '\r') ++i;
                    prevNonStandard = false;
	    	    break;
	    	case 'a': case 'A':
	    	    gSequence[gene_num++] = 'A';
                    prevNonStandard = false;
	    	    break;
	    	case 'c': case 'C':
	    	    gSequence[gene_num++] = 'C';
                    prevNonStandard = false;
	    	    break;
	    	case 'g': case 'G':
	    	    gSequence[gene_num++] = 'G';
                    prevNonStandard = false;
	    	    break;
	    	case 't': case 'T':
	    	    gSequence[gene_num++] = 'T';
                    prevNonStandard = false;
	    	    break;
	    	case '\n': case '\r':
	    	    ++blanks;
	    	    break;
	    	default: //Non-Standard nucleotides
                    if(!prevNonStandard){
	    	        gSequence[gene_num] = 'A'; 
	    	        long ranger = (gene_num + OLIG_LENGTH <= tempSize)? gene_num + OLIG_LENGTH : tempSize;
	    	        for (long j=gene_num; j< ranger; ++j)
	    	        	simarray[j] = 1;
	    	        ++errorno;
	    	        ++gene_num;
                        prevNonStandard = true;
                    } else{
                        curNcor++;
                        if(curNcor == 1 || NcorPos[curNcorIndex - 1] != gene_num+1){
                            if(curNcorIndex >= NcorSize){
                                NcorSize *= 2;
                                NcorPos = (long*) realloc(NcorPos, NcorSize * sizeof(long));
                                Ncor = (int*) realloc(Ncor, NcorSize * sizeof(int));
                            }
                            NcorPos[curNcorIndex] = gene_num+1;
                            Ncor[curNcorIndex] = curNcor;
                            curNcorIndex++;
                        } else{
                            Ncor[curNcorIndex-1]++;
                        }
                    }
	    }
	    ++i;
	}

	free(tempSequence);
	memoryUsed -= (tempSize) * sizeof(unsigned char);

	gSize = gene_num;
	if (gSize < tempSize){
	    gSequence = (unsigned char *) realloc(gSequence, gSize * sizeof(unsigned char));
	    simarray = (uint8_t *) realloc(simarray, gSize * sizeof(uint8_t));	
	}
        
	memoryUsed -= tempSize * (sizeof(unsigned char) + sizeof(uint8_t));
	memoryUsed += gSize * (sizeof(unsigned char) + sizeof(uint8_t));//resize gSeq and simarray
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	
	if (geneno+1 < locSize){
            geneLoc = (long*) realloc(geneLoc, (geneno+1) * sizeof(long));
            geneLen = (int *) realloc(geneLen, geneno * sizeof(int));
            geneNcor = (int *) realloc(geneNcor, geneno * sizeof(int));
        }
        geneLoc[geneno] = gSize;
	memoryUsed += (geneno+1) * sizeof(long);//geneLoc
	memoryUsed += (2*geneno) * sizeof(int);//geneLen + geneNcor

        if(curNcorIndex <= NcorSize){
            NcorPos = (long*) realloc(NcorPos,curNcorIndex * sizeof(long));
            Ncor = (int*) realloc(Ncor,curNcorIndex * sizeof(int));
        }
        NcorSize = curNcorIndex;
        memoryUsed += NcorSize * sizeof(long);//NcorPos
        memoryUsed += NcorSize * sizeof(int);//Ncor
	
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;

        //printf("NcorSize:%d\n",NcorSize);
        /*printf("***\n");
        printf("NcorPos\tNcor\n");
        for(int i = 0; i < NcorSize; ++i){
            printf("%ld\t%d\n",NcorPos[i],Ncor[i]);
        }
        printf("***\nGene\tNcor\n");
        for(int i = 0; i < geneno; i++){
            printf("%d\t%d\n",i,geneNcor[i]);
        }*/
	
	gseqOff = gSize % 4;
	
	if(gseqOff != 0){
		for(int i=1; i < geneno+1; ++i)
			geneLoc[i] -= gseqOff;
		for(long i=0; i < gSize - gseqOff; ++i)
			simarray[i] = simarray[i+gseqOff];
                geneLen[0] -= gseqOff;
	}
	
	nbytes = gSize/4;
        ////////////// Encode The Genome
        gCoded = (uint8_t *) malloc (gSize/4 * sizeof(uint8_t));
        if(!encodeSequence(gCoded,gSize,gSequence)) {
            cout << "Could not encode blacklist sequence\n";
            exit(EXIT_FAILURE);
        }
	memoryUsed += nbytes * sizeof(uint8_t);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
    	gSize -= gseqOff;
	simarray = (uint8_t *) realloc(simarray, gSize * sizeof(uint8_t));
        memoryUsed -= gseqOff * sizeof(uint8_t);
        for (int j = 0; j < geneno; j++){
	    for (long i = geneLoc[j], k = 0; k < OLIG_LENGTH - 1; ++i, ++k){
	        simarray[i] = 1;
                //simarray[i + geneLen[j]] = 1;
            }
        }
	//////////////////////////////// END create genome	
        
        printf("\nGenome size is %ldbp from %d target sequences\n",(gSize + gseqOff)/*/2*/, geneno);
        printf("%d Non-standard nucleotides encoded as Adenine/Thymine pairs\n",errorno);
        print_Count(simarray,geneLoc);
        //print_simArray(simarray);

	unsigned long memoryEval = 0;
        
	// Kane Identity Check 	
	PRIME_MM = pickPrime();
	memoryUsed += PRIME_MM * (sizeof(uint64_t) + sizeof(long));//hashTable and posTable
        memoryUsed += gSize * (sizeof(long) + sizeof(long*));//association count and table
	
	timers = time(NULL);	
	memoryEval = identity_Check(gCoded, simarray, geneLoc,nbytes);
        timerf = time(NULL);
        memoryUsed += (memoryEval * sizeof(long)); //number of associated positions
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
        printf("Identity checked in %ld sec\n", timerf - timers);
        printf("Step used %luMB of memory \n",((memoryUsed/1024)/1024));
        print_Count(simarray,geneLoc);
        //print_simArray(simarray);
	
	memoryUsed -= PRIME_MM * (sizeof(uint64_t) + sizeof(long));//hashTable and posTable
        memoryUsed -= memoryEval * sizeof(long); //associations
        memoryUsed -= gSize * (sizeof(long) + sizeof(long*));//association count and table
	
 	MAX_NUM_THREAD = omp_get_max_threads();	
        if(maxThread != -1 and maxThread < MAX_NUM_THREAD){
            MAX_NUM_THREAD = maxThread;
        }

	// GC Content Check
        timers = time(NULL);
	gc_Content(gCoded, simarray,nbytes);
        timerf = time(NULL);
        printf ("GC Content check done in %ld sec\n",(timerf-timers));
        print_Count(simarray,geneLoc);
	
	// Tm Check
	timers = time(NULL);
	float *tmarray = (float *) malloc (gSize*sizeof(float));
	memoryUsed += (gSize) * sizeof(float);
	if (memoryUsed > peakUsed) peakUsed = memoryUsed;
	for (long i = 0; i < gSize; ++i) 
		tmarray[i] = 0;
	melting_Temp(gCoded, simarray, geneLoc, tmarray);
	timerf = time(NULL);
        printf("Peak memory used in MB: %lu \n",((memoryUsed/1024)/1024));
        print_Count(simarray,geneLoc);	

	
        //Intensive Homology Search
	timers = time(NULL);
        memoryEval = 0;
        OLIGO* OLIGList = NULL;
	memoryEval = intensive_Homology(gCoded, simarray, gSequence, geneLoc, tmarray, geneLen,NcorSize,NcorPos,Ncor,geneNcor,OLIGList);
	if (memoryUsed  + memoryEval > peakUsed) peakUsed = memoryUsed + memoryEval;
        memoryUsed += sizeof(OLIGO) * olig_total;
	timerf = time(NULL);
	printf ("Intensive homology search done! in %ld sec\n", (timerf-timers));
        printf("Step used %luMB of memory\n",((memoryEval/1024)/1024));
        //print_Count(simarray,geneLoc);
        printf("There are %d candidate oligos from %d gene\n\n",olig_total,geneno);
        //print_simArray(simarray);

	free(gCoded);
        free(tmarray);
	free(simarray);	
	free(geneLoc);
        free(geneLen);
        free(geneNcor);
        free(NcorPos);
        free(Ncor);

        /////Extract each filename from the given list of blacklist files and process
        if(blackliststr != NULL){
            char* commaPtr = NULL;
            int commaPos = 0;

	    timers = time(NULL);
            while((commaPtr = strchr(&blackliststr[commaPos],',')) || commaPos < strlen(blackliststr)){
                int nameLen = ((commaPtr == NULL) ? strchr(blackliststr,'\0') : commaPtr) - &blackliststr[commaPos];
                unsigned char* file;

                file = (unsigned char*) malloc(sizeof(unsigned char) * nameLen + 1); 
                substr(file,(unsigned char*)(&blackliststr[commaPos]),strlen(blackliststr),0,nameLen);
                commaPos += nameLen + 1;

                memoryEval = ProcessBlacklistFile((char*)file,OLIGList);
	        if (memoryUsed  + memoryEval > peakUsed) peakUsed = memoryUsed + memoryEval;
                free(file);
            }
	    timerf = time(NULL);
	    printf ("Blacklisting Done! in %ld sec\n", (timerf-timers));
            printf("Step used used %luMB of memory\n\n",((memoryEval/1024)/1024));
        }

	oligos = fopen(oligosname, "w");
        for(int i = 0; i < olig_total; ++i){
            print_Oligo(OLIGList[i],oligos);
        }
	fclose(oligos);
	
	secondf = time (NULL);
        printf("Total number of unique oligos are: %d for %d genes.\n",olig_total, geneno);
        printf("Peak memory used in MB: %lu \n",((peakUsed/1024)/1024));
	printf ("Running Time: %ld sec\n\n", (secondf-seconds));	
	//for(i=0;i<olig_total;i++){
        //    free(OLIGList[i].sequence);
        //}

	free(gSequence);
        free(OLIGList);
	cout<<endl;
	print_Time();
	cout << "------------------------------END--------------------------------------\n";
    return 0;
}
