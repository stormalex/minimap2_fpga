#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

struct ref_name_t {
    char name[256];
    int  rid;
};

struct ref_name_t ref_name[1237514];
int    sort_name[1237514];

static int cmpstringp(const void *p1, const void *p2){
    return strcmp(((struct ref_name_t *)p1)->name, ((struct ref_name_t *)p2)->name);
}

unsigned int dichotomy_sort(char *qname, struct ref_name_t* ref_name, int ref_name_size)
{
    int i;
    unsigned int start=0, end = ref_name_size-1, mid; 
    int cmp;
    while(start < end) 
    {    
        mid = (start + end)>>1;
        cmp = strcmp(qname, ref_name[mid].name);
        if(cmp == 0)
            return (mid | 1UL<<31);
        else if( cmp < 0 )
            end = mid;
        else
            start = mid+1;
    }   
    return start;
}

int main()
{
    FILE *fp;
    char buff[256];
    int read_count = 0;
    int tmp;
    int i,j;
    clock_t t_start, t_end;

    //char test_search[]=">m150314_042626_42199_c100751742550000001823157807081544_s1_p0/108976/14848_19246";
    //char test_search[]=">m150314_042626_42199_c100751742550000001823157807081544_s1_p0/108976/14848_19246\n";
    char test_search[]="@8f7ff4cb-aeac-49e3-ac2a-41d82be67566_Basecall_Alignment_template PLSP57501_20161122_FNFAB45321_MN16458_sequencing_run_Hu_Bir_R94_1Dlig_fc2_1hr_81526_ch18_read1720_strand\n";
    unsigned int return_pos;
    
    //read reference name
    //fp = fopen("read_name2","r");
    fp = fopen("read_name1","r");
    if(fp == NULL)
    {
        printf("Open file eror!\n");
        return 0;
    }
    memset(buff, 0, 256);
    t_start = clock();
    while( fgets(buff, 256, fp) != NULL)
    {
        memcpy(ref_name[i].name, buff, 256);
        ref_name[i].rid = read_count;
        i++;
        read_count++;
        memset(buff, 0, 256);
    }
    fclose(fp);
    t_end = clock();
    printf("read file time %f\n", (t_end - t_start + 0.0)/CLOCKS_PER_SEC);
    fflush(stdout);

    //create table1
    t_start = clock();
    qsort(ref_name,read_count,sizeof(struct ref_name_t),cmpstringp);
    t_end = clock();
    printf("create table1 time %f\n", (t_end - t_start + 0.0)/CLOCKS_PER_SEC);
    fflush(stdout);

    //create table2
    t_start = clock();
    for(i=0;i<read_count;i++)
        sort_name[ref_name[i].rid] = i;
    t_end = clock();
    printf("create table2 time %f\n", (t_end - t_start + 0.0)/CLOCKS_PER_SEC);
    fflush(stdout);


    //display
    printf("\n\ntable1\n");
    for(i=0;i<read_count;i++)
        printf("%06d %s",ref_name[i].rid, ref_name[i].name);

#if 0
    printf("\n\ntable2\n");
    for(i=0;i<read_count;i++)
        printf("%08d\n", sort_name[i]);
    fflush(stdout);
#endif
    return_pos = dichotomy_sort(test_search, ref_name, read_count);
    printf("++++++++++++ %x %x %d %s\n", return_pos, return_pos>>31, return_pos & 0x7fffffff, test_search);
    fflush(stdout);
    return 0;
}
