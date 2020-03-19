#!/usr/bin/perl 
use strict;
use Getopt::Long;
use vars qw($snpiden $fasta $db $output $cv %subjects @segid %seen %count  @matsub %missub  @perfect @imperfect  @intrap @intraimp @interp @interimp);

GetOptions(
   "h|help"=>\&USAGE,
   "fsi=s"=>\$snpiden,
   "fa=s"=>\$fasta,
   "db:s"=>\$db,
   "o=s"=>\$output,
   "cv=s"=>\$cv,
) or USAGE();

sub USAGE{
my $usage=<<EOF;
Usage:perl $0  -fsi flank_snp.identity.txt  -fa wtf.fasta -db dbname  -o output  -cv copy_number.var.txt 
-h|--help:    print manual;
-fsi:   file recording pair-wise identies of flanking SNPs, which should involve 4 columns:strain1,strain2,identity(upstream),identity(downstream)
Format example:
JB938	JB22	0.99	1
Notice: when strainA vs strainB has occurred,strainB vs strainA is accepeted,but not recommended;	
-fa:    file containing all wtf sequnces in all strains in one wtf locus;
-db:    database_name for blast to search donors;
-o:	output
-cv:	copy number varations
EOF

print $usage;
exit;
};

`:>$output.tmp`;
open IDEN,$snpiden;
open CV,">$cv";
while(<IDEN>){
   chomp;
   my @arr=split;
   if($arr[2]<0.95 || $arr[3]<0.95){ next; }
   my$strain1 = $arr[0];
   my$strain2 = $arr[1];
   my$n1=`cat $fasta|grep $strain1|wc -l`;
   my$n2=`cat $fasta|grep $strain2|wc -l`;
   if($n1 == 0 || $n2==0){ next;}
   if($n1 != $n2){ 
      print CV "Recent copy number varation for $strain1 and $strain2\n";
      next;
   }
   my@name=split(/\n/,`cat $fasta|grep $strain1|sed 's/>//'|awk '{split(\$1,a,"_"); print a[2];}'`);
   for my $iname(@name){
       push(@{$subjects{"${strain1}_${iname}"}},"${strain2}_${iname}");
       push(@{$subjects{"${strain2}_${iname}"}},"${strain1}_${iname}");
       my$seq1=uc(`cat  $fasta|awk 'BEGIN{RS=">";OFS="";}{if(\$1=="${strain1}_${iname}"){\$1="";printf \$0;}}'`);
       my$seq2=uc(`cat  $fasta|awk 'BEGIN{RS=">";OFS="";}{if(\$1=="${strain2}_${iname}"){\$1="";printf \$0;}}'`);
       if($seq1 eq $seq2){ next;}

       `echo -e ">${strain1}_${iname}\n$seq1\n>${strain2}_${iname}\n$seq2" >${strain1}_vs_${strain2}.$iname.fa`;
       `mafft --quiet ${strain1}_vs_${strain2}.$iname.fa >${strain1}_vs_${strain2}.$iname.mafft.fa`;
       my$alnseq1=`cat  ${strain1}_vs_${strain2}.$iname.mafft.fa|awk 'BEGIN{RS=">";OFS="";}{if(\$1=="${strain1}_${iname}"){\$1="";printf \$0;}}'`;    
       my$alnseq2=`cat  ${strain1}_vs_${strain2}.$iname.mafft.fa|awk 'BEGIN{RS=">";OFS="";}{if(\$1=="${strain2}_${iname}"){\$1="";printf \$0;}}'`;
       my $posi1=0;
       my $posi2=0;
       open BED,">${strain1}_vs_${strain2}.$iname.mismatch.bed";
       for my $i(1..length($alnseq1)){
          my $nt1=substr($alnseq1,$i-1,1);
          if($nt1 ne "-"){ $posi1++; }else{};
          my $nt2=substr($alnseq2,$i-1,1);
          if($nt2 ne "-"){ $posi2++; }else{};
          if($nt1 ne "-" && $nt2 ne "-" && $nt1 ne $nt2){
              print BED "${strain1}_${iname}\t",$posi1-1,"\t$posi1\t${strain2}_${iname}\n${strain2}_${iname}\t",$posi2-1,"\t$posi2\t${strain1}_${iname}\n"; 
          }
       }
       close BED;
       if(-z "${strain1}_vs_${strain2}.$iname.mismatch.bed"){
       }else{
          `cat ${strain1}_vs_${strain2}.$iname.mismatch.bed | sort -k1,1 -k2,2n |  bedtools merge -d 198 -c 1,4 -o count,distinct -i -|awk 'BEGIN{OFS="\t"}{if(\$4>5 && (\$3-\$2)/\$4<=100){\$2++; print \$0;}}' >>$output.tmp`;#each 2 mismatches within 200bp are classified into one segment and the number of mismatches in a candidate GC segment should be >5;
       }
       `rm  ${strain1}_vs_${strain2}.$iname.fa ${strain1}_vs_${strain2}.$iname.mafft.fa ${strain1}_vs_${strain2}.$iname.mismatch.bed`;
   }
}
close IDEN;
close CV;

if(-z $cv){`rm $cv`;}
if(-z "$output.tmp"){ print "Done!No gene conversion detected.\n";`rm $output.tmp`; exit;}

###############duplicates remove and  largely overlap segments merging########
open TMP,"cat $output.tmp|sort -k1,1 -k2,2n -k3,3nr |";
open EXT,">$output.seg_extract.txt";
$_=<TMP>;
my @arr=split;
my $gene=$arr[0];
my $start=$arr[1];
my $end=$arr[2];
@{$missub{"$gene:$start-$end"}}=$arr[4];
while(<TMP>){
    chomp;
    @arr=split;
    if($arr[0] eq $gene){
        if($arr[2]<=$end){
             push(@{$missub{"$gene:$start-$end"}},$arr[4]);
        }elsif($end-$arr[1]+1>=0.8*($end-$start+1)){
             @{$missub{"$gene:$start-$arr[2]"}}=(@{$missub{"$gene:$start-$end"}},$arr[4]); 
             $end=$arr[2];            
        }else{
             print EXT "$gene\t$start\t$end\t+\t$gene:$start-$end\n";
             $start=$arr[1];
             $end=$arr[2];
             @{$missub{"$gene:$start-$end"}}=$arr[4];
        }
    }else{
        print EXT "$gene\t$start\t$end\t+\t$gene:$start-$end\n";
        $gene=$arr[0];
        $start=$arr[1];
        $end=$arr[2];
        @{$missub{"$gene:$start-$end"}}=$arr[4];
    }
}
print EXT "$gene\t$start\t$end\t+\t$gene:$start-$end\n";
close TMP;
close EXT;

###########similar segments filtering################
`perl /data/XYH/script/03_extract_fasta_from_ref.pl  $fasta  $output.seg_extract.txt  $output.seg.fa`;
@segid=split(/\n/,` cat $output.seg_extract.txt|awk '{print \$5;}'`);

for my$i(0..$#segid-1){
    if($seen{$segid[$i]}==1){next;}
    for my$j($i+1..$#segid){
        if($seen{$segid[$j]}==1){next;}
        $segid[$i]=~/(\S+):/; my$name1=$1;
        $segid[$j]=~/(\S+):/; my$name2=$1;
        if($name1 eq $name2){ next;}
        my$seg1=uc(`cat  $output.seg.fa|awk 'BEGIN{RS=">";OFS="";}{if(\$1=="$segid[$i]"){\$1="";printf \$0;}}'`); 
        my$seg2=uc(`cat  $output.seg.fa|awk 'BEGIN{RS=">";OFS="";}{if(\$1=="$segid[$j]"){\$1="";printf \$0;}}'`);
        if($seg1 eq $seg2){ $seen{$segid[$j]}=1;next;}
        `echo -e ">$segid[$i]\n$seg1\n>$segid[$j]\n$seg2" >$segid[$i]_vs_$segid[$j].fa`;
        `mafft --quiet $segid[$i]_vs_$segid[$j].fa >$segid[$i]_vs_$segid[$j].mafft.fa`;
        my $alnseg1=`cat $segid[$i]_vs_$segid[$j].mafft.fa|awk 'BEGIN{RS=">";OFS="";}{if(\$1=="$segid[$i]"){\$1=="";printf \$0;}}'`;
        my $alnseg2=`cat $segid[$i]_vs_$segid[$j].mafft.fa|awk 'BEGIN{RS=">";OFS="";}{if(\$1=="$segid[$j]"){\$1=="";printf \$0;}}'`;
        my $gaphead1; my $gaptail1; my $gaphead2; my $gaptail2;
        $alnseg1=~/^(-*)\w+\S+\w+(-*)$/;$gaphead1=$1;$gaptail1=$2;
        $alnseg2=~/^(-*)\w+\S+\w+(-*)$/;$gaphead2=$1;$gaptail2=$2;
        if(length($gaphead1)>=0.2*length($seg2) || length($gaptail1)>=0.2*length($seg2) || length($gaphead2)>=0.2*length($seg1) || length($gaptail2)>=0.2*length($seg1)){ `rm $segid[$i]_vs_$segid[$j].fa $segid[$i]_vs_$segid[$j].mafft.fa`; next;}
        my $mismatch=0;
        for my $i(1..length($alnseg1)){
            my $nt1=substr($alnseg1,$i-1,1);
            my $nt2=substr($alnseg2,$i-1,1);
            if($nt1 ne "-" && $nt2 ne "-" && $nt1 ne $nt2){ $mismatch++;}
        }
        if( $mismatch<=5 || (length($seg1)-length($gaphead2)-length($gaptail2))/$mismatch >100 || (length($seg2)-length($gaphead1)-length($gaptail1))/$mismatch >100 ){$seen{$segid[$j]}=1;}
        `rm $segid[$i]_vs_$segid[$j].fa $segid[$i]_vs_$segid[$j].mafft.fa`; 
    }
}

#`rm $output.seg_extract.txt $output.seg.fa`;

#######################produce final output#######
open SEG,"<$output.seg_extract.txt";
open OUT,">$output";
print OUT "gene\tGC_start\tGC_end\t#consistent_relatives\tconsistent_relatives\t#inconsistent_relatives\tinconsistent_relatives\tperfect_intradonors\timperfect_intradonors\tperfect_interdonors\timperfect_interdonors\n";
while(<SEG>){
    chomp;
    my @arr=split;
    if($seen{"$arr[0]:$arr[1]-$arr[2]"}==1){ next;}
##########consistent relatives#######
    $arr[0]=~/(\S+)_(\S+)/;
    my $istrain=$1;
    my $iname=$2;
    `:>$arr[0].subjects.fa`;
    ` cat  $fasta|awk 'BEGIN{RS=">";}{if(\$1=="$arr[0]")print ">"\$0;}'| grep -v '^\$' >>$arr[0].subjects.fa`;
    for my$isub(@{$subjects{$arr[0]}}){
       ` cat  $fasta|awk 'BEGIN{RS=">";}{if(\$1=="$isub")print ">"\$0;}'| grep -v '^\$' >>$arr[0].subjects.fa`;
    }
    `mafft --quiet  $arr[0].subjects.fa >$arr[0].subjects.mafft.fa`;
    `perl /data/XYH/script/fasta_algnment_cut.pl $arr[0].subjects.mafft.fa  $arr[0]:$arr[1]-$arr[2]|sed 's/\\/[0-9]\\+-[0-9]\\+//' >$arr[0].subjects.cut.fa`;
    my $iseq=` cat $arr[0].subjects.cut.fa|awk 'BEGIN{RS=">";}{if(\$1~"$arr[0]") printf \$2;}' `;
    @matsub=split(/\n/,`cat $arr[0].subjects.cut.fa|awk 'BEGIN{RS=">";}{if(\$1!~"$istrain" && \$2=="$iseq") print \$1;}'`);

########dornor search######
    @intrap=();@interp=();@intraimp=();@interimp=();
    $iseq=~s/-//g;
    ` echo -e  ">$arr[0]\n$iseq" >$arr[0]_$arr[1]-$arr[2].query.fa`;
    ` blastn -query $arr[0]_$arr[1]-$arr[2].query.fa  -db $db -outfmt 7  -out query-$arr[0]_$arr[1]-$arr[2].out`;
    @perfect=split(/\n/,`cat query-$arr[0]_$arr[1]-$arr[2].out|awk '{split(\$1,a,"_");split(\$2,b,"_");if(a[2]!=b[2] && \$3==100 && \$4==$arr[2]-$arr[1]+1 ) print \$2;}'`);
    for my$iperf(@perfect){
        $iperf=~/(\S+)_(\S+)/;
        if($iperf=~/$istrain/ || grep{ $_ =~ $1 }@matsub){ push(@intrap,$iperf);}else{push(@interp,$iperf);}
    }
    #intrap:perfect intra-genome donors; interp:perfect inter-genome donors;
    @imperfect=split(/\n/,`cat query-$arr[0]_$arr[1]-$arr[2].out|awk '{split(\$1,a,"_");split(\$2,b,"_");if(a[2]!=b[2] && \$3>99 && \$4>0.99*($arr[2]-$arr[1]+1) ) print \$2;}'`);
     %count=();@imperfect = grep {$count{$_} == 1 } grep {++$count{$_}} (@perfect, @imperfect);
     for my$iimper(@imperfect){
         $iimper=~/(\S+)_(\S+)/;
         if($iimper=~/$istrain/ || grep{ $_ =~ $1 }@matsub){ push(@intraimp,$iimper);}else{push(@interimp,$iimper);}
     }
     #intraimp:imperfect intra-genome donors; interimp:imperfect inter-genome donors;
     %count=();@{$missub{"$arr[0]:$arr[1]-$arr[2]"}}=grep { ++$count{ $_ } < 2; }@{$missub{"$arr[0]:$arr[1]-$arr[2]"}};
     print OUT join("\t",@arr[0..2]),"\t",$#matsub+1,"\t",join(",",@matsub),"\t",$#{$missub{"$arr[0]:$arr[1]-$arr[2]"}}+1,"\t",join(",",@{$missub{"$arr[0]:$arr[1]-$arr[2]"}}),"\t",join(",",@intrap),"\t",join(",",@intraimp),"\t",join(",",@interp),"\t",join(",",@interimp),"\n";
    `rm $arr[0].subjects.fa $arr[0].subjects.mafft.fa $arr[0].subjects.cut.fa $arr[0]_$arr[1]-$arr[2].query.fa query-$arr[0]_$arr[1]-$arr[2].out`;
}
close SEG;
close OUT;

#`rm $output.tmp`;
