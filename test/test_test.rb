#!/usr/bin/env	ruby

require 'helper'

class TestRecipSnpFinder < Test::Unit::TestCase

  context "RecipSnpFinder" do

    setup do
      # this is run before each test
      @l1 = BetterSam.new("FCC00CKABXX:2:1101:1599:2070#CAGATCAT	83	nivara_3s	1979736	42	100M	=	1979623	-213	ACGATGAGATGGCTAGGGATTGGGCACGAGAAGAGGAACGCCGCGGCGAGGCGGACGGTCGAGAGCTTCGTCGCGTCGGCCATCGCGAAACACAGAGCCG	BBBBBBBBBBbb^`T``][[SY[ZGUKYH\_Tdae__Z\`]ae_deVddggbgeceggeegegagegeggeg_gfggggggeefgcggdggggggggggg	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:100	YS:i:-23	YT:Z:CP")
      @l2 = BetterSam.new("FCC00CKABXX:2:1101:1599:2070#CAGATCAT	163	nivara_3s	1979623	42	95M1I4M	=	1979736	213	TCGACCCCGGCTGCCTGGACATCGGCCTGCCGGAGGTCCCCTTCGCACGGGCGATGGACGACGTGCTGCGGACGATCTTCCCCCGCCACCCCATGCCCCG	gggggggggegeggggggggbefddgga^gceeZeaYYbbb[b]\V]YW][IQS^TXZ[Z`[_U[b^M\c^J^[M]]\[^BBBBBBBBBBBBBBBBBBBB	AS:i:-23	XN:i:0	XM:i:3	XO:i:1	XG:i:1	NM:i:4	MD:Z:81T7A7A1	YS:i:0	YT:Z:CP")

      @list_of_snps = []
      @list_of_snps << Vcf.new("nivara_3s	779750	779750	A	C	15	Q20;GOF	BRF=0.23;FR=0.5000;HP=1;HapScore=1;MGOF=52;MMLQ=33;NF=4;NR=0;PP=15;RMSMQ=39.52;SC=CCGCGGCAGCAGCATCATGAA;SbPval=1.0;Source=Platypus;TC=24;TCF=22;TCR=2;TR=4	GT:GL:GOF:GQ:NR:NV	0/1:-5.58,0.0,-84.18:52:56:24:4")
      @list_of_snps << Vcf.new("nivara_3s	879750	8979750	A	C	15	Q20;GOF	BRF=0.23;FR=0.5000;HP=1;HapScore=1;MGOF=52;MMLQ=33;NF=4;NR=0;PP=15;RMSMQ=39.52;SC=CCGCGGCAGCAGCATCATGAA;SbPval=1.0;Source=Platypus;TC=24;TCF=22;TCR=2;TR=4	GT:GL:GOF:GQ:NR:NV	0/1:-5.58,0.0,-84.18:52:56:24:4")
      @list_of_snps << Vcf.new("nivara_3s	979750	979750	A	C	15	Q20;GOF	BRF=0.23;FR=0.5000;HP=1;HapScore=1;MGOF=52;MMLQ=33;NF=4;NR=0;PP=15;RMSMQ=39.52;SC=CCGCGGCAGCAGCATCATGAA;SbPval=1.0;Source=Platypus;TC=24;TCF=22;TCR=2;TR=4	GT:GL:GOF:GQ:NR:NV	0/1:-5.58,0.0,-84.18:52:56:24:4")
      @list_of_snps << Vcf.new("nivara_3s	1979750	1979750	A	C	15	Q20;GOF	BRF=0.23;FR=0.5000;HP=1;HapScore=1;MGOF=52;MMLQ=33;NF=4;NR=0;PP=15;RMSMQ=39.52;SC=CCGCGGCAGCAGCATCATGAA;SbPval=1.0;Source=Platypus;TC=24;TCF=22;TCR=2;TR=4	GT:GL:GOF:GQ:NR:NV	0/1:-5.58,0.0,-84.18:52:56:24:4")
      @list_of_snps << Vcf.new("nivara_3s	2979750	2979750	A	C	15	Q20;GOF	BRF=0.23;FR=0.5000;HP=1;HapScore=1;MGOF=52;MMLQ=33;NF=4;NR=0;PP=15;RMSMQ=39.52;SC=CCGCGGCAGCAGCATCATGAA;SbPval=1.0;Source=Platypus;TC=24;TCF=22;TCR=2;TR=4	GT:GL:GOF:GQ:NR:NV	0/1:-5.58,0.0,-84.18:52:56:24:4")
      @list_of_snps << Vcf.new("nivara_3s	3979750	3979750	A	C	15	Q20;GOF	BRF=0.23;FR=0.5000;HP=1;HapScore=1;MGOF=52;MMLQ=33;NF=4;NR=0;PP=15;RMSMQ=39.52;SC=CCGCGGCAGCAGCATCATGAA;SbPval=1.0;Source=Platypus;TC=24;TCF=22;TCR=2;TR=4	GT:GL:GOF:GQ:NR:NV	0/1:-5.58,0.0,-84.18:52:56:24:4")

      #def initialize            (origin,         target,          left,          right,   cores,  output, test, verbose)
      @snpfinder = SnpFinder.new "sativa.fasta", "barthii.fasta", "pooled.l.fq", "pooled.r.fq", 4, "output", false, true
    end

    should "find a snp" do
      assert @l1.contains_snp?(1979737) == true
    end

    should "not find a snp" do
      assert @l1.contains_snp?(1979730) == false
    end

    should "do a binary search" do
      assert @l1.find_snp(@list_of_snps) == 3
    end

    should "find correct snp" do
      
    end
  end
end