#!/usr/bin/perl
#
# DTR2NetLogo.pl
#
# Routine to process the output of a decision tree and output it as a netlogo script
#
#
# Example rpartScore decision tree:
#
# n=105 (17 observations deleted due to missingness)
#
# node), split, n, deviance, yval
#       * denotes terminal node
#
#  1) root 105 160 6  
#    2) BH_MOBILITY< 6.5 95 113 6  
#      4) IdentityEO< 6.5 46  58 6  
#        8) BH_MOBILITY< 5.125 38  39 6  
#         16) VBiospheric< 4.625 12  18 5 *
#         17) VBiospheric>=4.625 26  17 6 *
#        9) BH_MOBILITY>=5.125 8   9 4 *
#      5) IdentityEO>=6.5 49  55 6  
#       10) VBiospheric< 6.375 29  17 7  
#         20) BW_ENERGY_CONSUMPTION< 6.225 19  13 6 *
#         21) BW_ENERGY_CONSUMPTION>=6.225 10   1 7 *
#       11) VBiospheric>=6.375 20  31 5  
#         22) BW_WASTE_MANAGEMENT< 6.0625 8  11 3 *
#         23) BW_WASTE_MANAGEMENT>=6.0625 12  14 7 *
#    3) BH_MOBILITY>=6.5 10   3 1 *
# 
#
# Example rpart decision tree with method = "class":
#
# n=487 (131 observations deleted due to missingness)
#
# node), split, n, loss, yval, (yprob)
#       * denotes terminal node
#
#   1) root 487 283 3 (0.15400411 0.28952772 0.41889117 0.13757700)  
#     2) IdentityES< 4.166667 91  37 2 (0.17582418 0.59340659 0.16483516 0.06593407) *
#     3) IdentityES>=4.166667 396 207 3 (0.14898990 0.21969697 0.47727273 0.15404040)  
#       6) NormsP< 5.375 160 104 3 (0.18750000 0.30000000 0.35000000 0.16250000)  
#        12) VHedonic>=2.166667 153 104 3 (0.19607843 0.31372549 0.32026144 0.16993464)  
#          24) as.factor(Ed_Level)=3,4 78  54 1 (0.30769231 0.21794872 0.30769231 0.16666667)  
#            48) Worldviews< 6.416667 53  32 1 (0.39622642 0.24528302 0.24528302 0.11320755)  
#              96) NormsP< 4.875 34  17 1 (0.50000000 0.26470588 0.14705882 0.08823529)  
#               192) ES>=4.166667 27  11 1 (0.59259259 0.14814815 0.18518519 0.07407407) *
#               193) ES< 4.166667 7   2 2 (0.14285714 0.71428571 0.00000000 0.14285714) *
#              97) NormsP>=4.875 19  11 3 (0.21052632 0.21052632 0.42105263 0.15789474) *
#            49) Worldviews>=6.416667 25  14 3 (0.12000000 0.16000000 0.44000000 0.28000000)  
#              98) NormsIG< 3.5 8   1 3 (0.00000000 0.12500000 0.87500000 0.00000000) *
#              99) NormsIG>=3.5 17  10 4 (0.17647059 0.17647059 0.23529412 0.41176471) *
#          25) as.factor(Ed_Level)=2,5 75  44 2 (0.08000000 0.41333333 0.33333333 0.17333333)  
#            50) NormsIL< 3.166667 18   8 3 (0.16666667 0.22222222 0.55555556 0.05555556) *
#            51) NormsIL>=3.166667 57  30 2 (0.05263158 0.47368421 0.26315789 0.21052632)  
#             102) IWO>=4.166667 29  11 2 (0.03448276 0.62068966 0.27586207 0.06896552) *
#             103) IWO< 4.166667 28  18 4 (0.07142857 0.32142857 0.25000000 0.35714286)  
#               206) EfficacySelf< 4.5 14   7 2 (0.00000000 0.50000000 0.35714286 0.14285714) *
#               207) EfficacySelf>=4.5 14   6 4 (0.14285714 0.14285714 0.14285714 0.57142857) *
#        13) VHedonic< 2.166667 7   0 3 (0.00000000 0.00000000 1.00000000 0.00000000) *
#       7) NormsP>=5.375 236 103 3 (0.12288136 0.16525424 0.56355932 0.14830508)  
#        14) VBiospheric< 8.875 226  95 3 (0.11946903 0.17256637 0.57964602 0.12831858) *
#        15) VBiospheric>=8.875 10   4 4 (0.20000000 0.00000000 0.20000000 0.60000000) *
#
# The program expects an output on which the previous line contains the name of a behaviour
#
# Usage: DTR2NetLogo.pl <R output file> <NetLogo output file>
#
# Version 2, 1 December 2013
# For discretised trees (names ending in 'disc'), use the entries in a 'discrete' file
# which is also output by buildTrees2.R, to sample in a range. This is only done if
# the -discrete option is given.
#
# Version 3, 1 December 2013
# Add an option to produce 'DOT' graphs
#
# Version 4, 3 December 2013
# Add an option to produce a netlogo '__include' .nls file
# 
# Version 5, 19 May 2014
# Change the way the discretisation adding/subtracting 1e-6 is corrected for
# as not all behaviours have a minimum of 1 or a maximum of 7. buildTrees4.R
# does not add or subtract 1e-6, so this is no longer necessary.
use strict;

# Globals

my $default_node_chr = ')';
my $default_terminal_chr = '*';
my $default_split_info = 'split';
my $default_yval_info = 'yval';
my $default_n_info = 'n';
my $default_root_str = 'root';
my @default_info = ($default_split_info, $default_n_info, 'deviance/loss', $default_yval_info);

my $level_label = '_level';
my $terminal_label = '_terminal';
my $nodeid_label = "_nodeid";

my $default_netlogo_var_func = 'get-var';
my $default_netlogo_value_func = 'get-value';

my %discrete;
my %trees;
my $dotdir;
my $nls = 0;

# Process command line

while($ARGV[0] =~ /^-/) {
  my $opt = shift(@ARGV);

  if($opt eq "-discrete") {
    my $file = shift(@ARGV);

    open(FP, "<", $file) or die "Cannot open discrete file $file: $!\n";

  LINE:
    while(my $line = <FP>) {
      $line =~ s/\s*$//;

      my @cells = split(/,/, $line);

      next if length($cells[0]) == 0 or $cells[0] =~ /^\#/;

      my $behav = shift(@cells);

      foreach my $cell (@cells) {
	next LINE if $cell eq "NA";
      }

      # Correct for what was done in the R script. Minimum and maximum value
      # must always be an integer. As of version 5, this is no longer done,
      # but the commented out code below is still a change from version 4.
      #
      # $cells[0] = sprintf("%.0f", $cells[0]);	
      # $cells[$#cells] = sprintf("%.0f", $cells[$#cells]);

      $discrete{"$behav.disc"} = \@cells;
      $discrete{"$behav.disc.union"} = \@cells;
      $discrete{"$behav.disc.all"} = \@cells;
    }

    close(FP);
  }
  elsif($opt eq "-dotdir") {
    $dotdir = shift(@ARGV);

    if(-e "$dotdir" && ! -d "$dotdir") {
      die "Dot directory $dotdir exists and is not a directory\n";
    }
    elsif(!-e "$dotdir") {
      mkdir $dotdir or die "Cannot create dot directory $dotdir: $!\n";
    }
  }
  else {
    die "Option $opt not recognised\n";
  }

}

if(scalar(@ARGV) != 2) {
  die "Usage: $0 <R output file> <NetLogo output file>\n";
}

my $Rinput = shift(@ARGV);
my $NLoutput = shift(@ARGV);

if($NLoutput =~ /\.nls$/) {
  $nls = 1;
}

open(RIN, "<", $Rinput) or die "Cannot open R output file $Rinput: $!\n";

open(NLOUT, ">", $NLoutput) or die "Cannot create NetLogo output file $NLoutput: $!\n";

if($nls) {
  print NLOUT "to-report build-trees\n  let tree-table table:make\n";
}

my $prev_line;
my $comment_info = "NA";
my %behav;
my $all_comment = 1;

while(my $line = <RIN>) {
  $line =~ s/\s+$//;

  if($line =~ /Tree for behaviour(.+)with formula/) {
    $comment_info = $1;
    $comment_info =~ s/\s//g;
  }

  if($line =~ /^n(\s*)=(\s*)(\d+)\s*(\(.+\))?/) {
    my $n = $3;
    my $missing = 0;
    my $ok = 1;
    if(defined($4)) {
      my $missingstr = $4;
      if($missingstr =~ /^\((\d+)/) {
	$missing = $1;
      }
      else {
	warn "Ignoring line \"$line\" as it doesn't look like the start of a tree\n";
	$ok = 0;
      }
    }

    if($ok) {
      my $tree_name = $prev_line;
      my $behav_name = $prev_line;
      my $tree_id = $prev_line;

      if(defined($trees{$tree_name})) {
	$trees{$tree_name} += 1;
      }
      else {
	$trees{$tree_name} = 1;
      }

      $behav_name =~ s/\.disc$//;

      if($comment_info ne "NA") {
	if($comment_info =~ /union/) {
	  $tree_id .= ".union";
	}
	elsif($comment_info =~ /all/) {
	  $tree_id .= ".all";
	}
      }
      else {
	$tree_id .= ".NA$trees{$tree_name}";
	$all_comment = 0;
      }

      push(@{$behav{$behav_name}}, $tree_id);

      print NLOUT "; tree found $prev_line $n observations $missing missing\n";

      if($nls) {
	print NLOUT "  table:put tree-table \"$tree_id\"";
      }
      else {
	print NLOUT "$tree_id\n";
      }

      my $tree = &read_tree(*RIN, $n, []);

      if(defined($dotdir)) {
	open(DOT, ">", "$dotdir/$tree_id.dot")
	  or die "Cannot create dot file $dotdir/$tree_id.dot: $!\n";
	&print_tree($tree_id, *NLOUT, $tree, *DOT);
	close(DOT);
      }
      else {
	&print_tree($tree_id, *NLOUT, $tree);
      }
    }
    $comment_info = "NA";
  }
  else {
    # Ignore the line
  }

  $prev_line = $line if $line !~ /^\s*$/;
}

if($nls) {
  print NLOUT "  report tree-table\nend\n";
}

close(RIN);
close(NLOUT);

if($all_comment) {
  print "behav,normal,discrete,normal.union,discrete.union,normal.all,discrete.all\n";

  foreach my $b (sort keys(%behav)) {
    print $b;
    
    my %which;
    
    for(my $i = 0; $i < scalar(@{$behav{$b}}); $i++) {
      $which{$behav{$b}->[$i]} = 1;
    }

    foreach my $suffix ("", ".disc", ".union", ".disc.union", ".all", ".disc.all") {
      if(defined($which{"$b$suffix"})) {
	print ",$b$suffix";
      }
      else {
	print ",NA";
      }
    }
    
    print "\n";
  }
}
else {
  foreach my $b (sort keys(%behav)) {
    print $b, ",", join(",", @{$behav{$b}}), "\n";
  }
}

exit 0;

#################################################################################################

sub read_tree {
  my ($fp, $n, $tree) = @_;

  my ($terminalchr, $nodechr, @header) = &read_tree_header($fp);

  my $nterminal = 0;		# Tree done when $nterminal >= $n (should be ==)
  while(my $line = <$fp>) {
    next if $line =~ /^\s*$/;	# Ignore empty lines

    $line =~ s/\s*$//;

    if($line =~ /^(\s*)(\d+)(.*)$/) {

      my $spaces = $1;
      my $nodeid = $2;
      my $rest = $3;

      if(substr($rest, 0, 1) ne $nodechr) {
	warn "Not found node character $nodechr where expected on line \"$line\" -- ignoring\n";
	next;
      }

      my $terminal = substr($rest, -1) eq $terminalchr;

      if($terminal) {
	$rest = substr($rest, 2, -2);
				# Expect spaces after nodechr and before terminalchr
      }
      else {
	$rest = substr($rest, 2);
				# Not a terminal node
      }
      
      my $var;
      my $op;
      my $value;

      if($rest =~ /[<>=]/) {	# There's a comparison operator in there somewhere...
	my @varopdata = split(/([<>= ])/, $rest);

	$var = shift(@varopdata);
				# Assume the first thing is the variable

	if($var =~ /^as\.factor\((.+)\)$/) {
	  $var = $1;		# Remove as.factor from node split
	}

	# Build the operator
	while($varopdata[0] =~ /[<>= ]/ || $varopdata[0] eq "") {
	  my $chr = shift(@varopdata);
	  $op .= $chr if $chr =~ /[<>=]/;
	}

	$value = shift(@varopdata);

	$rest = join("", @varopdata);
      }
      else {			# Should be a root node because there's no comparison
	if(substr($rest, 0, 4) eq $default_root_str) {
	  $rest = substr($rest, 4);
	  $rest =~ s/^\s*//;
	}
	else {
	  warn "Line \"$line\" ought to be a root node, but I don't detect it as such -- ",
	    "ignoring\n";
	  next;
	}

	$var = "root";
	$op = "";
	$value = "";
      }
    

      # By this point, $var, $op and $value should contain what's in the 'split' part of
      # the info, and $rest should contain everything else. We expect this everything
      # else to be separated by spaces, except where there's a bracket.
    
      my %node;			# This is the node data

      $node{$default_split_info} = [$var, $op, $value];
      $node{$terminal_label} = $terminal;
      $node{$level_label} = length($nodeid) + length($spaces);
      $node{$nodeid_label} = $nodeid;

      my @restsep = split("[() ]", $rest);

      for(my $i = 0; $i <= $#header; $i++) {
	next if $header[$i] eq $default_split_info;

	# Remove spaces and nulls
	while($restsep[0] =~ /^\s*$/) {
	  shift(@restsep);
	}

	if($restsep[0] eq "(") {
				# This one is an array
	  my @arr;

	  # Remove spaces and nulls
	  while($restsep[0] =~ /^\s*$/) {
	    shift(@restsep);
	  }

	  while($restsep[0] ne ")" && scalar(@restsep) > 0) {
	    push(@arr, shift(@restsep));

	    # Remove spaces and nulls
	    while($restsep[0] =~ /^\s*$/) {
	      shift(@restsep);
	    }

	  }

	  $node{$header[$i]} = \@arr;
	}
	else {			# Not an array
	  $node{$header[$i]} = shift(@restsep);
	}
      }

      $nterminal += $node{$default_n_info} if $terminal;

      push(@$tree, \%node);

      last if $nterminal >= $n;
    }
    else {
      # print "Ignoring line $line\n";
    }
  }

  if($nterminal != $n) {
    warn "Terminal nodes added up to $nterminal instead of $n\n";
  }

  return $tree;
}

sub read_tree_header {
  my ($fp) = @_;

  my @header;
  my $chr;
  my $nodechr;
  while(my $line = <$fp>) {
    next if $line =~ /^\s*$/;

    if($line =~ /^\s+(.)\s+denotes\s+terminal\s+node\s+$/i) {
      # The line is telling us about the terminal node indicator
      # We are flexible about this character, which is * at present in R
      $chr = $1;
      last;
    }
    else {
      # Assume the line is telling us metadata about each line in the tree
      # We are flexible about $nodechr, which is ) at present in R
      my @poss_headers = split(/,/, $line);

      for(my $i = 0; $i <= $#poss_headers; $i++) {
	$poss_headers[$i] =~ s/^\s*//;
	$poss_headers[$i] =~ s/\s*$//;
      }

      if($poss_headers[0] =~ /^node(.)$/i) {
        $nodechr = $1;
	shift(@poss_headers);
	@header = @poss_headers;
      }
    }
  }

  if(!defined($chr) || !defined($nodechr) || scalar(@header) == 0) {
    warn("Not found all metadata about the tree\n");
    $chr = $default_terminal_chr if !defined($chr);
    $nodechr = $default_node_chr if !defined($nodechr);
    @header = @default_info if scalar(@header) == 0;
  }

  return($chr, $nodechr, @header);
}

sub print_tree {
  my ($name, $fp, $tree, $dotfp) = @_;

  my %levels;
  for(my $i = 0; $i <= $#$tree; $i++) {
    push(@{$levels{$tree->[$i]->{$level_label}}}, $i);
				# For each level, build an array of
				# nodes as indexes of @$tree in the
				# order they appear
  }

  my @levs = sort { $a <=> $b } (keys(%levels));

  if(defined($dotfp)) {
    print $dotfp "digraph tree {\n  node \[fontname = \"Helvetica\", fontsize = 14.0, ",
      "shape = \"diamond\"\];\n  edge \[fontname = \"Helvetica\", fontsize = 12.0, ",
	"labelfontname = \"Helvetica\", labelfontsize = 10.0, arrowhead = \"vee\", ",
	  "arrowtail = \"vee\"\];\n";
  }

  &print_node($name, $fp, \@levs, \%levels, 0, $tree, 0, $dotfp);

  if(defined($dotfp)) {
    print $dotfp "}\n";
  }
}

sub print_node {
  my ($name, $fp, $levels, $level_nodes, $i, $tree, $nodeid, $dotfp) = @_;

  my $nodearr = $level_nodes->{$$levels[$i]};

  if($i == 0) {			# Root node, and there should be just one of them
    if(scalar(@$nodearr) > 1) {
      print STDERR "There seems to be more than one root node (";

      for(my $j = 0; $j <= $#$nodearr; $j++) {
	print STDERR ", " if $j > 0;
	print STDERR $tree->[$j]->{$nodeid_label};
      }

      print STDERR ") -- only the first will be used\n";
    }

    my $root = shift(@$nodearr);

    if($tree->[$root]->{$default_split_info}->[0] ne $default_root_str) {
      warn "What should be a root node doesn't have the $default_split_info variable name ",
	"\"$default_root_str\" (\"", $tree->[$root]->{$default_split_info}->[0], "\")\n";
    }

    if($nls) {
      print $fp " task [";
    }
    else {
      print $fp "report ";
    }

    &print_node($name, $fp, $levels, $level_nodes, $i + 1, $tree, $tree->[$root]->{$nodeid_label}, $dotfp);

    if($nls) {
      print $fp "]";
    }

    print $fp "\n";
  }
  else {			# Non-root node, and they should always be in pairs
    if(!defined($nodearr)) {
      warn "No defined nodes at level $i\n";
      print $fp "false";

      return;
    }

    if(scalar(@$nodearr) < 2) {
      warn "Expecting at least two nodes remaining at level $i (got ", scalar(@$nodearr), ")\n";
      print $fp "false";

      return;
    }

    my $lhs = shift(@$nodearr);
    my $rhs = shift(@$nodearr);

    print $fp "ifelse-value (";

    my $lhs_test = $tree->[$lhs]->{$default_split_info};

    my ($lhs_var, $lhs_op, $lhs_value) = @$lhs_test;

    my $rhs_test = $tree->[$rhs]->{$default_split_info};

    my ($rhs_var, $rhs_op, $rhs_value) = @$rhs_test;

    my $lhs_id = $tree->[$lhs]->{$nodeid_label};
    my $rhs_id = $tree->[$rhs]->{$nodeid_label};

    if($lhs_var ne $rhs_var) {
      warn "What should be a fork is comparing different variables \"$lhs_var\" and ",
	\"$rhs_var\" (nodes ", $tree->[$lhs]->{$nodeid_label}, " and ",
	  $tree->[$rhs]->{$nodeid_label}, ")\n";
    }

    if(defined($dotfp)) {
      my $lhsdv = ($lhs_value =~ /,/) ? "\[$lhs_value\]" : sprintf("%.3g", $lhs_value);
      my $rhsdv = ($rhs_value =~ /,/) ? "\[$rhs_value\]" : sprintf("%.3g", $rhs_value);
      print $dotfp "  node_$nodeid \[label = \"\", width = 0.25, height = 0.25, regular = true\];\n";
      print $dotfp "  node_$nodeid -> node_$lhs_id \[label = \"\[$lhs_var $lhs_op $lhsdv\]\"\];\n";
      print $dotfp "  node_$nodeid -> node_$rhs_id \[label = \"\[$rhs_var $rhs_op $rhsdv\]\"\];\n";
    }
    
    if($lhs_op ne '=' && $rhs_op ne '=' && $lhs_value ne $rhs_value) {
      warn "Inequality tests $lhs_op and $rhs_op are not comparing on the same values ",
	"($lhs_value and $rhs_value) at what should be a fork in the tree (nodes ",
	  $tree->[$lhs]->{$nodeid_label}, " and ", $tree->[$rhs]->{$nodeid_label}, ")\n";
    }

    if(($lhs_op eq ">" && $rhs_op ne "<=")
       || ($lhs_op eq ">=" && $rhs_op ne "<")
       || ($lhs_op eq "<" && $rhs_op ne ">=")
       || ($lhs_op eq "<=" && $rhs_op ne ">")) {
      warn "Tests $lhs_op and $rhs_op are not complementary at what should be a fork in ",
	"the tree (nodes ", $tree->[$lhs]->{$nodeid_label}, " and ",
	  $tree->[$rhs]->{$nodeid_label}, ")\n";
    }

    if($tree->[$lhs]->{$level_label} != $tree->[$rhs]->{$level_label}) {
      die "Fork at nodes with different levels (nodes ",
	$tree->[$lhs]->{$nodeid_label}, " with level ", $tree->[$lhs]->{$level_label},
	  " and ", $tree->[$rhs]->{$nodeid_label}, " with level ",
	    $tree->[$rhs]->{$level_label}, ")\n"
    }

    my @lhs_values = split(/,/, $lhs_value);
				# For =, there may be a list of (factor) values

    # Print the test

    for(my $j = 0; $j <= $#lhs_values; $j++) {
      print $fp " or " if $j > 0;
      print $fp "(" if scalar(@lhs_values) > 1;
      print $fp "$default_netlogo_var_func \"$lhs_var\" $lhs_op ",
	"$default_netlogo_value_func \"$lhs_var\" $lhs_values[$j]";
      print $fp ")" if scalar(@lhs_values) > 1;
    }

    print $fp ") [";

    # Print the value if the test is true

    if($tree->[$lhs]->{$terminal_label}) {
      my $value = $tree->[$lhs]->{$default_yval_info};
      if(defined($discrete{$name})) {
	my $min = $discrete{$name}->[$value - 1];
	my $max = $discrete{$name}->[$value];

	my $mind = sprintf("%.3g", $min);
	my $maxd = sprintf("%.3g", $max);

	if($max == $min) {
	  print $fp $max;

	  if(defined($dotfp)) {
	    print $dotfp "  node_$lhs_id \[shape = \"rect\", style = \"rounded\", ",
	      "label = \"$maxd\"\];\n";
	  }
	}
	else {
	  print $fp "$min + random-float ", ($max - $min);

	  if(defined($dotfp)) {
	    print $dotfp "  node_$lhs_id \[shape = \"rect\", style = \"rounded\", ",
	      "label = \"\[$mind, $maxd\]\"\];\n";
	  }
	}
      }
      else {
	print $fp $value;

	if(defined($dotfp)) {
	  print $dotfp "  node_$lhs_id \[shape = \"rect\", style = \"rounded\", ",
	    "label = \"", sprintf("%.3g", $value), "\"\];\n";
	}
      }
    }
    else {
      print $fp " ";
      &print_node($name, $fp, $levels, $level_nodes, $i + 1, $tree, $lhs_id, $dotfp);
      print $fp " ";
    }

    print $fp "] [";

    # Print the value if the test is false

    if($tree->[$rhs]->{$terminal_label}) {
      my $value = $tree->[$rhs]->{$default_yval_info};
      if(defined($discrete{$name})) {
	my $min = $discrete{$name}->[$value - 1];
	my $max = $discrete{$name}->[$value];

	my $mind = sprintf("%.3g", $min);
	my $maxd = sprintf("%.3g", $max);

	if($max == $min) {
	  print $fp $max;
	  if(defined($dotfp)) {
	    print $dotfp "  node_$rhs_id \[shape = \"rect\", style = \"rounded\", ",
	      "label = \"$maxd\"\];\n";
	  }
	}
	else {
	  print $fp "$min + random-float ", ($max - $min);
	  if(defined($dotfp)) {
	    print $dotfp "  node_$rhs_id \[shape = \"rect\", style = \"rounded\", ",
	      "label = \"\[$mind, $maxd\]\"\];\n";
	  }
	}
      }
      else {
	print $fp $value;
	if(defined($dotfp)) {
	  print $dotfp "  node_$rhs_id \[shape = \"rect\", style = \"rounded\", ",
	    "label = \"", sprintf("%.3g", $value), "\"\];\n";
	}
      }
    }
    else {
      print $fp " ";
      &print_node($name, $fp, $levels, $level_nodes, $i + 1, $tree, $rhs_id, $dotfp);
      print $fp " ";
    }

    print $fp "]";
  }
}
