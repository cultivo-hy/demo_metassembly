#!/usr/bin/env python
import os
f1 = open('../data/SRR800798_1.fastq', 'r')
f2 = open('../data/SRR800798_1', 'w')
bool = 1;
for line in f1:
    if(bool % 2 == 1) :
        line = line[:11] + str(bool/4 +1) +'/1' + '\r\n'
    f2.write(line)
    bool +=1
print 'finish'
f1.close()
f2.close()

os.system('perl -p -i -e "s/\r//g" /home/ubuntu/data/SRR800798_1')

f1 = open('../data/SRR800798_2.fastq', 'r')
f2 = open('../data/SRR800798_2', 'w')
bool = 1;
for line in f1:
    if(bool % 2 == 1) :
        line = line[:11] + str(bool/4 +1) +'/2' + '\r\n'
    f2.write(line)
    bool +=1
print 'finish'
f1.close()
f2.close()

os.system('perl -p -i -e "s/\r//g" /home/ubuntu/data/SRR800798_2')
