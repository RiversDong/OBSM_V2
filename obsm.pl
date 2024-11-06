#!/usr/bin/perl -w
use open ':std', ':utf8';
use Cwd;
use File::Spec;
use Getopt::Long;
use File::Path qw(make_path);
use File::Copy;
use Cwd 'abs_path';
use File::Basename;
my $main_dir = dirname(abs_path(__FILE__));
require "$main_dir/utils.pl";

=pop
Version: 2.0
Authors: Chuan Dong, Ruiyuan Li and Chengjun Zhang#
Email: zhangcj@zafu.edu.cn

Please feel free to contact us with any comments or questions.
=cut


my $output_path = './result';
my $aln_file;
my $tree_file;
GetOptions(
    'o=s'     => \$output_path,
    'aln=s'   => \$aln_file,
    'tree=s'  => \$tree_file,
    'fo=s' => \$focus,
) or die "Cannot find the files!\n";

unless ($aln_file && $tree_file) {
    die "Error: You must specify the -aln and -tree parameters";
}
if ($output_path) {
    unless (-d $output_path) {
        make_path($output_path) or die "cannot creat $output_path: $!\n";
    }
}
if ($aln_file && !-e $aln_file) {
    die "there is no: $aln_file\n";
}
if ($tree_file && !-e $tree_file) {
    die "there is no: $tree_file\n";
}

unless ($focus && $focus ne '') {
    die "Error: The -fo (focus) parameter cannot be empty.\n";
}

my $codeml = File::Spec->catfile($main_dir, 'codeml');
#print $codeml."\n";
#print $main_dir."\n";
my $processor="16";
open (PERCENT,">percent.log");



&paml_analysis_main($output_path, $codeml);
close PERCENT;

=pop
  $path就是-o传来的参数
  &check_original_tree 是在-o 指定的路径下读取文件
 下面的$codeml 是codeml的绝对路径，是perl程序
 &get_ORM_rel($codeml);
 &get_FRM_rel("F", $codeml); 这两个程序其实可以并行 加快计算
 free ratio model get_FRM_rel 这个函数最为浪费时间最为浪费时间
=cut

sub paml_analysis_main {
	my ($path, $codeml)=@_;
	chdir $path;
	my $des_seq = File::Spec->catfile($path, 'dna_seq_for_paml.txt');
	copy($aln_file, $des_seq) or die "Failed to copy $aln_file to $des_seq: $!";
	my $des_tree = File::Spec->catfile($path, 'gene_tree.trees');
	copy($tree_file,  $des_tree) or die "Failed to copy $tree_file to $des_tree: $!";
	my $original_tree=&check_original_tree;
	&get_the_ctl_file("0","gene_tree.trees");
	&get_ORM_rel($codeml);
	#&get_the_ctl_file("1","gene_tree.trees");
	#&get_FRM_rel("F", $codeml); 
	# 物种和逗号分割成不同的行 $edit_tree
	# 共分割成了多少行
	# 摒弃了原来的程序停止的条件判断，改为主进程等待子进程执行结束 才接着往下运行
	&get_the_ctl_file("2","trm_tree.trees");
	my ($edit_tree, $edit_tree_lines)=&get_edit_tree($original_tree);
	my $labled_edit_tree;
	# 手动标注w不一样的分支，哪个标记表示不一样的分支？
	my $kkkkk=0;
	my @pids;
	for(my $i=0;$i<$edit_tree_lines;$i++) {
		$labled_edit_tree=&lable_one_anchor($edit_tree,$i);
		&get_trm_tree($labled_edit_tree);
		my $pid = &get_sbTRM_rel($i, $codeml);
		push @pids, $pid;
		$kkkkk=$i;
	}
	# {中间的进行了修改 等待所有的程序执行完成
	&get_the_ctl_file("1","gene_tree.trees");
	my $pid = &get_FRM_rel_v2("F", $codeml);
	push @pids, $pid;
	foreach my $pid (@pids) {
		waitpid($pid, 0);
	}
	print "I. Finish the calculation of one ration model\n";
	print "II. Finish the calculation of free ration model\n";
	print "III. Finish the calculation of multi-ratio model\n";
	# 中间的进行了修改}

	print "IV. Begin select the best model\n";
	# %表示定义了一个名为lnL_hash的哈希 
	# 下面的代码开始有问题了
	# 在 Perl 中，函数的返回值在列表上下文中将被视为一个数组，在标量上下文中则会返回数组的大小（即元素的个数）。如果在某些情况下，current_document_files 函数的返回值被解释为标量，你可能会看到这样的行为。
	# 如果你在调用函数时使用 & 符号（即 &current_document_files()），这会将其置于标量上下文中，可能导致返回值不符合预期。为了确保以列表上下文调用函数，请直接调用它,不要加上&
	my @second_level_ducoment = current_document_files();
	my $current_directory = `pwd`;
	#exit;
	
	# -o more_rm_times 为指定的子目录下面产生的文件的数量
	# 文件夹和文件夹下面的
	# lnL_hash 也有点问题
	# 在 Perl 中，如果你在调用函数时使用了 & 符号，返回的结果可能会被解释为标量，而不是预期的哈希。
	# 去掉 & 符号,调用函数时，不要使用 & 符号，这样可以确保以列表上下文调用函数，返回值将会被正确处理。
	my ($more_rm_times,%lnL_hash) = get_lnL_hash(@second_level_ducoment);
	$more_rm_times=$more_rm_times-4; #？
	my @new_pids;
	for (my $j=0;$j<$more_rm_times;$j++) {
		my $more_rm_tree=$edit_tree;
		my %temp_lnL_hash=%lnL_hash;
		my $more_rm_maker;
		for (my $k=-2; $k<$j; $k++){
			my @max_lnL=("-100000","");
			foreach my $key(sort keys %temp_lnL_hash){
				if (($key ne "FRM") and ($key ne "ORM") and ($max_lnL[0]<$temp_lnL_hash{$key})){
					$max_lnL[0]=$temp_lnL_hash{$key};
					$max_lnL[1]=$key;
					$max_lnL[2]=$temp_lnL_hash{$key}{'df'};
				}
			}
			delete ($temp_lnL_hash{$max_lnL[1]});			
			my $anchor_marker=$k+3;
			$more_rm_maker=$anchor_marker+1;
			$anchor_marker="#"."$anchor_marker";
			my ($temp,$anchor_line)=split/-/,$max_lnL[1];
			$more_rm_tree=&lable_more_anchor($more_rm_tree,$anchor_line,$anchor_marker);
		}
		&get_more_rm_tree($more_rm_tree);
		&get_the_ctl_file("2","more_rm_tree.trees");
		my $pid = &get_more_rm_rel($more_rm_maker, $codeml);
		push @new_pids, $pid;
		$kkkkk=$more_rm_maker;
	}
	foreach my $pid (@new_pids) {
		waitpid($pid, 0);
	}
	#print "finish get_more_rm_rel"."\n";
	#exit;

	# { 两者之间不删除rst是否可以 ？？
	for ($i=3;$i<=$kkkkk;$i++){
		while(&check_more_rm_finished($i)){}
		system "rm $i-RM/rst";
	}
	# }

	@second_level_ducoment=&current_document_files();
	&get_lnL_hash_with_document(@second_level_ducoment);
	system "rm -rf codeml.ctl more_rm_tree.trees trm_tree.trees";
	chdir "..";
}

my $evolution_force = "python $main_dir/best_model.py -p $output_path -fo $focus";
print $evolution_force;
system($evolution_force) == 0 or die "Failed to execute command: $!\n";