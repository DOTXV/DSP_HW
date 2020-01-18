#include "hmm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct{
    int first;
    double second;
}pair;

void init(pair* p){
    p->first = 0;
    p->second = 0.0;
}

double viterbi(const HMM* hmm, const char* seq){
    double de[MAX_SEQ][MAX_STATE]={{0.0}};
    int len = strlen(seq);

    for(int i=0;i<hmm->state_num;i++)
        de[0][i] = hmm->initial[i] * hmm->observation[seq[0]-'A'][i];

    for(int t=1;t<len;t++){
        for(int i=0;i<hmm->state_num;i++){
            double max = 0.0;
            for(int j=0;j<hmm->state_num;j++){
                double cur = de[t-1][j] * hmm->transition[j][i];
                if(cur > max)
                    max = cur;
            }
            de[t][i] = max * hmm->observation[seq[t]-'A'][i];
        }
    }

    double maxP = -1.0;
    for(int i=0;i<hmm->state_num;i++)
        if(maxP<de[len-1][i])
            maxP = de[len-1][i];
    
    return maxP;
}

void go(HMM* model, const char* seq, int n, pair* p){
    int bestM = -1;
    double maxP = 0.0;
    for(int i=0;i<n;i++){
        double curP = viterbi(model+i, seq);
        if(curP > maxP){
            maxP = curP;
            bestM = i;
        }
    
    }
    p->first = bestM;
    p->second = maxP;
}

int main(int argc, char const *argv[])
{
    HMM models[MAX_LINE];

    int model_cnt = load_models(argv[1], models, MAX_LINE);
    char seq[MAX_SEQ] ="";
    FILE* data = open_or_die(argv[2], "r");
    FILE* result = open_or_die(argv[3], "w");

    while(fscanf(data, "%s", seq)!=EOF){
        pair p;
        init(&p);
        go(models, seq, model_cnt, &p);
        //printf("%d %f", p.first, p.second);
        fprintf(result, "%s %e\n", models[p.first].model_name, p.second);
    }

    fclose(data);
    fclose(result);
    return 0;
}
