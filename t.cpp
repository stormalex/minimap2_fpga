#include <stdio.h>


#define RESIZE64(n) ((n+63)&(~63))

#define ALIGN_BYTE_N(n, v) ((v+n-1)&(~(n-1)))
int main(int argc, char *argv[])
{
for (int k=0; k<10; k++) printf("%d, 16B=%d, 64B=%d\n",k,ALIGN_BYTE_N(16,k), ALIGN_BYTE_N(64,k));
//return 0;

for (int x=0; x<50; x++) {
if (x%16==0) fprintf(stderr, "\nppp-ppp:");
fprintf(stderr, " %02x", x);
//if (x && x%16==0) fprintf(stderr, "\n");
}
fprintf(stderr, "\n");

return 0;
}
