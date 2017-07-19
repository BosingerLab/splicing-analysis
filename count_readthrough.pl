#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;

my %opts = ( 'verbose' => 1, 'extended' => 0, 'format' => "partek", 'readGroups' => "read_groups.info", 'tracking' => "genes.fpkm_tracking", 'quant' => "FPKM", 'label' => "Sample" );

GetOptions(\%opts, "verbose", "extended", "splicesites=s", "genes=s", "format=s", "annotation=s", "readGroups=s", "tracking=s", "xcise", "quant=s", "label=s");

my %results;

if($opts{splicesites}) {
    open(SS, "< $opts{splicesites}") or warn "Couldn't open splicesites file: $opts{splicesites}\n";
    while(<SS>) {
	chomp();
	my @fields = split(/\t/);
	my ($start, $end) = ($fields[0] =~ /(\d+)\<\-\>(\d+)/);
	if($start && $end) {
	    $results{splicelist}{$start}{$end} = $fields[1];
	}
    }
}

print "ID\t$opts{label}\n";
while(<>) {
    my @fields = split(/\t/);
    my $id = $fields[0];
    my $pos = $fields[3];
    my $cigar = $fields[5];
    if($opts{verbose} > 1) {
	print "$pos\t$cigar\n";
    }
#    $_ = $cigar;
    my $cur = "";
    my $il;
    while($cigar =~ /\G(\d+\w)/gc) {
	if($opts{verbose} > 1) {
	    print "$1\n";
	}
	my $comp = $1;
	my ($count,$type) = ($comp =~ /(\d+)(\w)/);
	if($opts{verbose} > 1) {
	    print "$count\t$type\n";
	}
	# if(($type ne "N") && ($type ne "M")) {
	#     print "$type\t$cigar\n";
	# }
	if($type eq "S") { next; }
	elsif($type eq "M") {
	    if($cur eq "intron") {
		my $end = $pos + $il;
#		print "$pos\t$end\n";
		$results{splices}{$pos}{$end}++;
	    } else {
## readthrough calculation
		if($results{splicelist}) {
		    $rtend = $pos + ($count -1);
		    for my $start (sort {$a <=> $b} keys %{$results{splicelist}}) {
			if(($start > $pos) && ($start < $rtend)) {
			    $results{readthrough}{$start}{count}++;
			    $results{readthrough}{$start}{list}{$id}++
			}
			for my $end (sort {$a <=> $b} keys %{$results{splicelist}{$start}}) {
			    if(($end > $pos) && ($end < $rtend)) {
				$results{readthrough}{$end}{count}++;
				$results{readthrough}{$end}{list}{$id}++
			    }
			}
		    }
		}
##
		$pos += ($count - 1);
	    }
	    $cur = "match";
	    $il = 0;
	} elsif($type eq "N") {
	    if($cur eq "match") {
		$cur = "intron";
		$il = $count;
	    }
	} elsif($type eq "I") {
	    $cur = "insertion";
	    $il = 0;
	} elsif($type eq "D") {
	    $cur = "deletion";
	    $il = 0;
	} else {
	    warn "Unhandled type $type\n";
	}
    }
}
if($results{splicelist}) {
    for my $start (sort {$a <=> $b} keys %{$results{splicelist}}) {
	for my $end (sort {$a <=> $b} keys %{$results{splicelist}{$start}}) {
	    if($results{splices}{$start}{$end}) {
		print "$start<->$end\t$results{splices}{$start}{$end}\n";
	    } else {
		print "$start<->$end\t0\n";
	    }
	}
    }
    my @readthroughlist;# = keys %{$results{splices}};
    for my $sloc (sort {$a <=> $b} keys %{$results{splicelist}}) {
	if(!grep(/$sloc/,@readthroughlist)) {
	    push(@readthroughlist,$sloc);
	}
	for my $eloc (sort {$a <=> $b} keys %{$results{splicelist}{$sloc}}) {
	    if(!grep(/$eloc/,@readthroughlist)) {
		push(@readthroughlist,$eloc);
	    }
	}
    }
#    for my $loc (sort {$a <=> $b} keys %{$results{readthrough}}) {
    for my $loc (sort {$a <=> $b} @readthroughlist) {
	if($results{readthrough}{$loc}) {
	    my $rc = keys %{$results{readthrough}{$loc}{list}};
#	    print "$loc\t$results{readthrough}{$loc}{count}\t$rc\n";
	    print "$loc\t$rc\n";
	} else {
#	    print "$loc\t0\t0\n";
	    print "$loc\t0\n";
	}
    }
} else {
    for my $start (sort {$a <=> $b} keys %{$results{splices}}) {
	for my $end (sort {$a <=> $b} keys %{$results{splices}{$start}}) {
	    print "$start<->$end\t$results{splices}{$start}{$end}\n";
	}
    }
    for my $loc (sort {$a <=> $b} keys %{$results{readthrough}}) {
	if($results{readthrough}{$loc}) {
	    my $rc = keys %{$results{readthrough}{$loc}{list}};
#	    print "$loc\t$results{readthrough}{$loc}{count}\t$rc\n";
	    print "$loc\t$rc\n";
	} else {
#	    print "$loc\t0\t0\n";
	    print "$loc\t0\n";
	}
    }
}
