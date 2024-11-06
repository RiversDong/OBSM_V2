#!/usr/bin/perl -w
use utf8;

sub check_original_file
{  
	if (!(-e "gene_tree.trees") or !(-e "dna_seq_for_paml.txt")){
		exit;
	}
}
1;

sub check_original_tree
{
	open(TREE,"<gene_tree.trees") or die "can't open the tree file";
	my $tree="";
	while(<TREE>)
	{
		$tree=$tree.$_;
	}
	close TREE;
	return $tree;
}

=pop
	model = 0 的时候   one ratio model
	model = 1 的时候   free ratio model
	model = 2 的时候   multi ratio model 需要手动标注w值不一样的分支
	这个函数主要是生成需要的condem的
=cut
sub get_the_ctl_file
{
	my ($model,$tree_file)=@_;
	open (CTL,">codeml.ctl");
	print CTL <<END;
	    seqfile = dna_seq_for_paml.txt
END
	
     print CTL "treefile = ",$tree_file,"\n";
  print CTL <<END;
      outfile = REL.txt
        noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

*       ndata = 10
        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = dat/mtArt.dat  * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
END
       print CTL "model = ",$model;
        print CTL <<END;
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0
                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1 * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 8  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
* fix_blength = -1 * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0  * Optimization method 0: simultaneous; 1: one branch a time

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
END
	
	close CTL;
	
}

=pop
 这个文件其实就是在输出路径下面再建立新的路径 然后运行codeml
 shift是获得第一个参数
 标准输入 (stdin)：文件描述符 0
 标准输出 (stdout)：文件描述符 1
 
=cut
sub get_ORM_rel {
	$codeml = shift;
	mkdir "ORM";
	system "cp gene_tree.trees dna_seq_for_paml.txt codeml.ctl ORM/";
	chdir "ORM";
	system("$codeml > /dev/null &") == 0 or die "Failed to execute $codeml: $!";
	chdir "..";
}

sub get_FRM_rel_v2 {
    my ($com_twice, $codeml, @more_lnL) = @_;
    if ($com_twice eq "F") {
        mkdir "FRM" unless -d "FRM";  # 仅在目录不存在时创建
        system "cp gene_tree.trees dna_seq_for_paml.txt codeml.ctl FRM/";
        chdir "FRM";

        my $pid = fork();  # 创建子进程
        if (not defined $pid) {
            die "Failed to fork: $!";
        } elsif ($pid == 0) {
            # 子进程执行的命令
            exec("$codeml > /dev/null") or die "Failed to exec $codeml: $!";
        }
        chdir "..";  # 返回父进程的目录
        return $pid;  # 返回子进程的 PID
    }
}

sub get_FRM_rel {

	my ($com_twice, $codeml, @more_lnL)=@_;
	if ($com_twice eq "F"){
		mkdir "FRM";
		system "cp gene_tree.trees dna_seq_for_paml.txt codeml.ctl FRM/";
		chdir "FRM";
		system("$codeml > /dev/null &") == 0 or die "Failed to execute $codeml: $!";
		chdir "..";
	}
	#elsif($com_twice eq "S"){
	#	system "/data/chuand/fusion_gene/obsm/OBSM/codeml";
	#	my $first_rel=`grep lnL FRM/REL.txt`;
	#	my $second_rel=`grep lnL REL.txt`;
	#	if ($first_rel ne $second_rel){
	#		&get_FRM_rel("M",$first_rel,$second_rel);
	#	}else{
	#		system "rm 2NG.dN 2NG.dS 2NG.t 4fold.nuc rub lnf REL.txt rst rst1";
	#	}
	#}else{
	#	mkdir "FRM-TEMP";
	#		system "mv 2NG.dN 2NG.dS 2NG.t 4fold.nuc rub lnf REL.txt rst rst1 FRM-TEMP/";
	#	system "/data/chuand/fusion_gene/obsm/OBSM/codeml";
	#	my $third_rel=`grep lnL REL.txt`;
	#	my ($first_df,$first_lnL)=&get_lnL_value($more_lnL[0]);
	#	my ($second_df,$second_lnL)=&get_lnL_value($more_lnL[1]);
	#	my ($third_df,$third_lnL)=&get_lnL_value($third_rel);
	#	if ($first_lnL>$second_lnL){
	#		if ($third_lnL>=$first_lnL){
	#			system "rm 2NG.dN 2NG.dS 2NG.t 4fold.nuc rub lnf REL.txt rst rst1";
	#			system "rm -rf FRM-TEMP";
	#		}elsif ($third_lnL<=$second_lnL){
	#			system "rm 2NG.dN 2NG.dS 2NG.t 4fold.nuc rub lnf REL.txt rst rst1";
	#			system "rm -rf FRM";
	#			system "mv FRM-TEMP FRM";
	#		}else{
	#			system "mv 2NG.dN 2NG.dS 2NG.t 4fold.nuc rub lnf REL.txt rst rst1 FRM/";
	#			system "rm -rf FRM-TEMP";
	#		}
	#	}else{
	#		if ($third_lnL<=$first_lnL){
	#			system "rm 2NG.dN 2NG.dS 2NG.t 4fold.nuc rub lnf REL.txt rst rst1";
	#			system "rm -rf FRM-TEMP";
	#		}elsif ($third_lnL>=$second_lnL){
	#			system "rm 2NG.dN 2NG.dS 2NG.t 4fold.nuc rub lnf REL.txt rst rst1";
	#			system "rm -rf FRM";
	#			system "mv FRM-TEMP FRM";
	#		}else{
	#			system "mv 2NG.dN 2NG.dS 2NG.t 4fold.nuc rub lnf REL.txt rst rst1 FRM/";
	#			system "rm -rf FRM-TEMP";
	#		}
	#	}
	#}
}

sub get_edit_tree
{
	my ($tree)=@_;
	my ($i,$j);
	while($tree=~s/\s//){}
	while($tree=~s/\n//){}
	while($tree=~s/,/+\n/){$i++;}
	while($tree=~s/\)/=\n/){$j++;}
	while($tree=~s/\+/,/){}
	while($tree=~s/\=/)/){}
	$i=$i+$j;
	return ($tree,$i);
}

=pop
手动标注分支上不同的w值
=cut
sub lable_one_anchor {
	my ($tree,$line)=@_;
	my @tree=split/\n/,$tree;
	my $last_letter=substr($tree[$line],-1,1);
	if ($last_letter eq ","){
		$tree[$line]=~s/,/#1,/;
	}
	if ($last_letter eq ")") {
		$tree[$line]=~s/\)/#1\)/;
	}
	$tree=join("",@tree);
	return $tree;
}

=pop
">trm_tree.trees" 表示按照写的方式打开
=cut
sub get_trm_tree
{
	my ($tree)=@_;
	open (TRM_TREE,">trm_tree.trees") or die "something is wrong at sub get_trm_tree";
	print TRM_TREE $tree;
	close TRM_TREE;
}

sub get_sbTRM_rel {
    my ($i, $codeml) = @_;
    mkdir "TRM-$i";
    system("cp trm_tree.trees dna_seq_for_paml.txt codeml.ctl TRM-$i/");
    chdir "TRM-$i";

    my $pid = fork();  # 创建子进程
    if (not defined $pid) {
        die "Failed to fork: $!";
    } elsif ($pid == 0) {
        # 子进程执行的命令
        exec("$codeml > /dev/null") or die "Failed to exec $codeml: $!";
    }
    
    chdir "..";  # 父进程执行
    return $pid;  # 返回子进程的 PID
}


=pop
 通过检查对应文件夹下面REL.txt文件是不是存在
=cut
sub check_part_finsh
{
	my ($i)=@_;
	if (-e "TRM-$i/REL.txt"){
		my $temp_lnL=`grep "Time used" TRM-$i/REL.txt`;
		my $temp_lnL_value=`grep "lnL" TRM-$i/REL.txt`;
		my $temp_lnL_Dstree=`grep "TreeView" TRM-$i/REL.txt`;
		if ($temp_lnL eq "" and ($temp_lnL_value eq "" and $temp_lnL_Dstree eq "")){
			sleep (2);
			return 1;
		}else{
			return 0;
		}
	}else{
		return 1;
	}
}


sub check_trm_finished
{
	my ($i)=@_;
	if (-e "TRM-$i/REL.txt"){
		my $temp_lnL=`grep "Time used" TRM-$i/REL.txt`;
		my $temp_lnL_value=`grep "lnL" TRM-$i/REL.txt`;
		my $temp_lnL_Dstree=`grep "TreeView" TRM-$i/REL.txt`;
		if ($temp_lnL eq "" and ($temp_lnL_value eq "" and $temp_lnL_Dstree eq "")){
			sleep (20);
			return 1;
		}else{
			return 0;
		}
	}else{
		return 1;
	}
}

sub check_FRM_finished
{
	my ($i)=@_;
	if (-e "FRM/REL.txt"){
		my $temp_lnL=`grep "Time used" FRM/REL.txt`;
		my $temp_lnL_value=`grep "lnL" FRM/REL.txt`;
		my $temp_lnL_Dstree=`grep "TreeView" FRM/REL.txt`;
		if ($temp_lnL eq "" and ($temp_lnL_value eq "" and $temp_lnL_Dstree eq "")){
			sleep (2);
			return 1;
		}else{
			return 0;
		}
	}else{
		return 1;
	}
}

# -o指定的路径中的所有文件
sub current_document_files{
	my @current_document_files=<*>;
	return @current_document_files;
}

sub get_lnL_hash {
	my @second_level_document=@_;
	my %lnL_hash;
	my $documents=0;
	foreach (@second_level_document){
		next if $_=~/\./;
		$documents++;
		my $temp_lnL=`grep lnL $_/REL.txt`;
		($lnL_hash{$_}{'df'},$lnL_hash{$_})=&get_lnL_value($temp_lnL);#前者是df,后者是lnL
	}
	return ($documents,%lnL_hash);
}

sub lable_more_anchor {
	my ($tree,$line,$marker)=@_;
	my @tree=split/\n/,$tree;
	my $last_letter=substr($tree[$line],-1,1);
	if ($last_letter eq ","){
		$tree[$line]=~s/,/$marker,/;
	}
	if ($last_letter eq ")"){
		$tree[$line]=~s/\)/$marker\)/;
	}
	$tree=join("\n",@tree);
	return $tree;
}

sub get_more_rm_tree
{
	my ($tree)=@_;
	open (TRM_TREE,">more_rm_tree.trees")	or die "something is wrong at sub get_more_rm_tree";
	print TRM_TREE $tree;
	close TRM_TREE;
}

sub get_more_rm_rel {
    my ($i, $codeml) = @_;
    mkdir "$i-RM" unless -d "$i-RM";  # 仅在目录不存在时创建
    system "cp more_rm_tree.trees dna_seq_for_paml.txt codeml.ctl $i-RM/";
    chdir "$i-RM";

    my $pid = fork();  # 创建子进程
    if (not defined $pid) {
        die "Failed to fork: $!";
    } elsif ($pid == 0) {
        # 子进程执行的命令
        exec("$codeml > /dev/null") or die "Failed to exec codeml: $!";
    }

    chdir "..";  # 返回父进程的目录
    return $pid;  # 返回子进程的 PID
}

sub check_moreRM_part_finsh
{
	my ($i)=@_;
	if (-e "$i-RM/REL.txt"){
		my $temp_lnL=`grep "Time used" $i-RM/REL.txt`;
		my $temp_lnL_value=`grep "lnL" $i-RM/REL.txt`;
		my $temp_lnL_Dstree=`grep "TreeView" $i-RM/REL.txt`;
		if ($temp_lnL eq "" and ($temp_lnL_value eq "" and $temp_lnL_Dstree eq "")){
			sleep (2);
			return 1;
		}else{
			return 0;
		}
	}else{
		return 1;
	}
}

sub check_more_rm_finished
{
	my ($i)=@_;
	if (-e "$i-RM/REL.txt"){
		my $temp_lnL=`grep "Time used" $i-RM/REL.txt`;
		my $temp_lnL_value=`grep "lnL" $i-RM/REL.txt`;
		my $temp_lnL_Dstree=`grep "TreeView" $i-RM/REL.txt`;
		if ($temp_lnL eq "" and ($temp_lnL_value eq "" and $temp_lnL_Dstree eq "")){
			sleep (2);
			return 1;
		}else{
			return 0;
		}
	}else{
		return 1;
	}
}

sub get_lnL_hash_with_document
{
	my @second_level_document=@_;
	my %lnL_hash;
	open (LAST_REL,">last_rel.txt");
	print LAST_REL "Model\tdf\tlnL\n";
	foreach (@second_level_document){
		next if $_=~/\./;
		my $temp_lnL=`grep lnL $_/REL.txt`;
		my @abcd=&get_lnL_value($temp_lnL);
		print LAST_REL	"$_\t$abcd[0]\t$abcd[1]\n";
	}
	close LAST_REL;
}

sub get_lnL_value
{
	my ($temp_lnL_value)=@_;
	while ($temp_lnL_value=~s/  / /){}
	my @temp=split/ /,$temp_lnL_value;
	$temp[3]=~s/\)://;
	my @abc=@temp[3,4];
	return @abc;
}



sub current_document_path
{
	my $current_document_path=getcwd(); # 当前的工作路径
	return $current_document_path;
}





sub CDFChi2                                                        
{
	my ($x,$v)=@_;
	return (&CDFGamma($x,($v)/2.0,0.5));
}
	
sub CDFGamma
{
	my ($x,$alpha,$beta)=@_;
	return (&IncompleteGamma(($beta)*($x),$alpha,&LnGammaFunction($alpha)));
}
sub IncompleteGamma 
{
   my ($x, $alpha, $ln_gamma_alpha)=@_;
   my $i;
   my $p=$alpha;
   my $g=$ln_gamma_alpha;
   my $accurate=1e-8; 
   my $overflow=1e30;
   my $factor; 
   my $gin=0; 
   my $rn=0; 
   my $a=0;
   my $b=0;
   my $an=0;
   my $dif=0;
   my $term=0; 
   #$pn[6];

   return (0) if ($x==0) ;
   return (-1) if ($x<0 or $p<=0);

   $factor=exp($p*log($x)-$x-$g);   
   goto l30 if ($x>1 and $x>=$p) ;
#   /* (1) series expansion */
   $gin=1;  $term=1;  $rn=$p;
 l20:
   $rn++;
   $term*=$x/$rn;   
   $gin+=$term;
	 goto l20 if ($term > $accurate) ;
   $gin*=$factor/$p;
   goto l50;
 l30:
#   /* (2) continued fraction */
   $a=1-$p;   
   $b=$a+$x+1;
   $term=0;
   my @pn;
   $pn[0]=1;  
   $pn[1]=$x;  
   $pn[2]=$x+1;  
   $pn[3]=$x*$b;
   $gin=$pn[2]/$pn[3];
 l32:
   $a++;  
   $b+=2;  
   $term++;   
   $an=$a*$term;
   for ($i=0; $i<2; $i++) {
   	$pn[$i+4]=$b*$pn[$i+2]-$an*$pn[$i];
   }
   goto l35 if ($pn[5] == 0) ;
   $rn=$pn[4]/$pn[5];   
   $dif=abs($gin-$rn);
   goto l34 if ($dif>$accurate);
   goto l42 if ($dif<=$accurate*$rn) ;
 l34:
   $gin=$rn;
 l35:
   for ($i=0; $i<4; $i++){
   		$pn[$i]=$pn[$i+2];
   }
   goto l32 if (abs($pn[4]) < $overflow);
   for ($i=0; $i<4; $i++){
   	$pn[$i]/=$overflow;
  }
   goto l32;
 l42:
   $gin=1-$factor*$gin;

 l50:
   return ($gin);
}
	
sub LnGammaFunction 
{
	 my ($alpha)=@_;
   my $x=$alpha; 
   my $f=0; 
   my $z;

   if ($x<7) {
       $f=1;  
       $z=$x-1;
       while(++$z<7){
         $f*=$z;
       }
       $x=$z;   
       $f=-log($f);
   }
   $z = 1/($x*$x);
   return  $f + ($x-0.5)*log($x) - $x + .918938533204673+ (((-.000595238095238*$z+.000793650793651)*$z-.002777777777778)*$z+.083333333333333)/$x;  
}