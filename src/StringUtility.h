#ifndef __STRING_UTILITY_H__
#define __STRING_UTILITY_H__
#include <string.h>
#include <string>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
int mystrcmp(const char* a, const char* b);
int mytrim(const char* input, char** result);
int mystrcmp(const char* a, const char* b){ 
    int i=0,j=0;
    while(isblank(a[i])){
        i++;
    }
    while(isblank(b[j])){
        j++;
    }
    //printf("%c\n",a[i]);
    //printf("%c\n",b[j]);
    while(a[i]!='\0'&&b[j]!='\0'){
        if(a[i]!=b[j]){
            //printf("%c != %c\n",a[i],b[j]);
            return 1;
        }
        i++;
        j++;
    }
    if((a[i]=='\0')&&(b[j]!='\0')){
        //printf("%c != %c\n",a[i],b[j]);
        while(b[j]!='\0'){
            if(!isblank(b[j])){
                return 1;
            }
            j++;
        }
    }
    if((b[j]=='\0')&&(a[i]!='\0')){
        //printf("%d'%c'!=%d'%c'\n",i,a[i],j,b[j]);
        while(a[i]!='\0'){
            if(!isblank(a[i])){
                return 1;
            }
            i++;
        }
    }
    return 0;
}
int mytrim(const char *input, char **result){
    int begin,end,length;
    length=strlen(input);
    //printf("length=%d\n",length);
    if(length==0){
        (*result)=(char*)malloc(1*sizeof(char));
        (*result)[0]='\0';
        return 1;
    }
    for(int ii=0; ii < length;++ii){
        if(!isspace(input[ii]) || ii==(length-1)){
            begin=ii;
            break;
        }
    }
    if(begin==length-1){
        //printf("'%s' is empty\n",input);
        (*result)=(char*)malloc(1*sizeof(char));
        (*result)[0]='\0';
        return 1;
    }
    for(int ii=length-1; ii>=begin; --ii){
       if(!isspace(input[ii])){
           end=ii;
           break;
       } 
    }
    //printf("'%s' begin=%d end=%d\n",input,begin,end);
    (*result)=(char*)malloc((end-begin+1)*sizeof(char));
    (*result)[end-begin+1]='\0';
    strncpy((*result),&input[begin],end-begin+1);
    return 0;
}
/*
int main(){
    std::string a="CA ";
    char *out;
    char *test;
    mytrim(a.c_str(),&out);
    printf("in='%s'\n",a.c_str());
    printf("out='%s'\n",out);
    return 0;
}*/
#endif
