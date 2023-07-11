#/usr/bin/perl
# Usage: perl bsmap2cx.pl <bsmap.result> > <CX_report>
# Description: Convert bsmap result to CX_report

open (IN,"$ARGV[0]") || die $!;
<IN>;
while (<IN>){
        chomp;
        s/^chr//g;
        my @line = split /\t/;
        my $tri = $line[3];
        my $con = "CHH";
        if ($line[2] eq "-"){
                $tri =~ tr/ATGC/TACG/;
                $tri = reverse($tri);
        }
        my @context =  split(//,$tri,3);
        $tri = $context[2];
        @context = split(//,$tri,3);
        if ($context[1] eq "G"){
                $con = "CG";
        }elsif ($context[2] eq "G"){
                $con = "CHG";
        }

        my @CX = ($line[0], $line[1], $line[2], $line[6], $line[7]-$line[6], $con, "CGT");
        print join("\t",@CX)."\n";
}
