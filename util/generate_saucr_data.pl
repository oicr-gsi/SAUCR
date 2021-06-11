use strict;
use warnings;
use JSON::PP;
use Data::Dumper;
use File::Basename;


my %opts = ();
GetOptions(
	"fpr|f:s"        		=> \$opts{"fpr"},
    "refresh|r:s"			=> \$opts{"refresh"},
	"prefix|p:s"     		=> \$opts{"prefix"},
);
validate_options(\%opts);


### generate a time stamp
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $ts = sprintf ( "%04d%02d%02d:%02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec); 

print STDERR "parsing current file provenance report and saving metatdata from each fingerprint\n";


### extract all fingerprintCollector records from FPR, along with metadata
my $rpt="/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz";

### fields from FPR
#1 datetime
#2 projects
#8 identity
#14 library
#18 geo info
#19 run
#27 pool
#47 file
#57 limskey
my @recs=`cat ../FPR.txt | grep fingerprintCollector | cut -f 1,2,8,14,18,19,23,27,33,37,47,57`;chomp @recs;
my @keys=qw/date project identity library info run platform pool wfid wfrunid file limskey/;


(open my $META,">",$opts{prefix} . "_metadata.txt") || die "unable to open metadata file for output";
my @fields=qw/limskey project identity external_name library tissue_origin tissue_type group_id library_source_template_type prep_kit run platform pool file wfid wfrunid sites sites_nocovered sites_covered sites_wildtype sites_variant/;
print $META join("\t",@fields) . "\n";

### collection of fin files for jaccard calcluations
my %fin;

### temp

@recs=@recs[0..40];


for my $rec(@recs){
	
	### store record information in a hash
	my @val=split /\t/,$rec;
	my %h;
	@h{@keys}=@val;
	
	### bypass if 1. no fin file or 2. no file
	next unless(-e $h{file});
	next unless($h{file}=~/fin$/);
	
	### map the info key/vlue pairs into the hash
	map{
		my($key,$val)=split /=/,$_;
		$key=~s/geo_//;
		$h{$key}=$val;
	} split /;/,$h{info};
	
	### modificatiosn to values in the hash
	$h{pool}=~s/pool_name=//;
	
	### collect stats from the fingerprint on the number of sites, covered sites, variant sites.  from col 5 of the fingerprint
	my @flags=`cat $h{file} | cut -f 5`;chomp @flags;
	shift @flags;
	$h{sites}=scalar @flags;
	$h{sites_notcovered}=grep /N/,@flags;$h{sites_notcovered}=$h{sites_notcovered} || 0;
	$h{sites_covered}=$h{sites}-$h{sites_notcovered};
	$h{sites_wildtype}=grep/M/,@flags;$h{sites_wildtype}=$h{sites_wildtype} || 0;
	$h{sites_variant}=$h{sites_covered}-$h{sites_wildtype};

    ### store the fingeprint, indexing by limskey and noting the workflowrunid
	#$fin{$h{limskey}}={wfid=>$h{wfrunid},fin=>$h{file}};
	$fin{$h{wfrunid}}={limskey=>$h{limskey},fin=>$h{file},status=>"new"};
	
	### print to the meta file
	my @vals=map{ $h{$_} || ""} @fields;
	print $META join("\t",@vals) . "\n";
}
#print Dumper(%fin);<STDIN>;

### check records in the previous metadata report to asses which wfrundis are no longer used, and which are new
if($opts{refresh}){
	my $meta=$opts{refresh} . "_metadata.txt";
	if(! -e $meta){
	print STDERR "refresh metadata $meta not found\n";
	
}
if(-e $meta){
	my @recs=`cat $meta`;chomp @recs;
	my @headers=split /\t/,shift @recs;
	map{
		my %h;
		@h{@headers}=split /\t/,$_;
		my $wfrunid=$h{wfrunid};
		if($fin{$wfrunid}){
			$fin{$wfrunid}{status}="inplace"
		}else{
			$fin{$wfrunid}{status}="remove";
		}
	}@recs;
}


### add records for any new wfids

#exit;
#print Dumper(%fin);
#<STDIN>;



print STDERR "starting jaccard score generation\n";
# module load sample-fingerprinting
#jaccard_coeff_pair
(open my $JACCARD,">","saucr_jaccard.latest.txt");

### carry over records from the last jaccard with valid wfrunids in the pair
if(-e $jacc){
	(open my $LASTJACC,"<",$jacc);
	for my $rec(<$LASTJACC>){
		chomp $rec;
		#print "$rec";<STDIN>;
		my @f=split /\t/,$rec;
		my($id1,$id2)=split /_/,pop @f;
		#print "$id1 $id2";<STDIN>;
		#print Dumper($fin{$id1});<STDIN>;
		#print Dumper($fin{$id2});<STDIN>;
		
		unless($fin{$id1}{status} eq "remove" || $fin{$id2}{status} eq "remove"){
			print $JACCARD "$rec\n";
		}
	}
}


print Dumper(%fin);<STDIN>;

### get any new fin files
my @new;
map{
	push(@new,$_) if($fin{$_}{status} eq "new");
	
}keys %fin;

print Dumper(@new);


my @wfrunids=sort {$fin{$a}{limskey} ne $fin{$a}{limskey}} keys %fin;
for my $i(0..$#wfrunids){
	my $wfrunid1=$wfrunids[$i];
	my $lk1=$fin{$wfrunid1}{limskey};
	my $fin1=$fin{$wfrunid1}{fin};
	my $status1=$fin{$wfrunid1}{status};
	
	next  unless(-e $fin1);
	for my $j($i..$#wfrunids){
		my $wfrunid2=$wfrunids[$i];
		my $lk2=$fin{$wfrunid2}{limskey};
		my $fin2=$fin{$wfrunid2}{fin};
		my $status2=$fin{$wfrunid2}{status};
		next unless(-e $fin2);

		next unless($status1 eq "new" || $status2 eq "new");
		print "$lk1 $lk2\n";
	
		#print "jaccard_coeff_pair -id1 $lk1 -fin1 $fin1 -id2 $lk2 -fin2 $fin2";exit;
		my $json=`perl ../scripts/jaccard_coeff.pair.pl -id1 $lk1 -fin1 $fin1 -id2 $lk2 -fin2 $fin2`;
		#print "$json";<STDIN>;
		my $h=decode_json $json;
		my $jaccard=$$h{jaccard_score} || 0;
		my $covered=$$h{pair}{covered} || 0;
		my $wfid_pair=$wfrunid1 . "_" . $wfrunid2;
		print $JACCARD "$lk1\t$lk2\t$jaccard\t$covered\t$wfid_pair\n";
	}
}





sub validate_options{
	my ($opts)=@_;
	
	usage("Help requested.") if($opts{help});
	if(! $opts{fpr}){
		$opts{fpr}="/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz";
	}
	if(! -e $opts{fpr}){
		usage("ERROR : file provenance not found at $opts{fpr}");
	}
	
	if(!$opts{prefix}){
		$opts{prefix}="saucr";
	}
	
	if(!$opts{refresh}){
		$opts{refreshMeta}=$opts{refresh} . "_metadata.txt";
		if(! -e $opts{refreshMeta}){
			usage("ERROR : meta data file for refresh not found, $opts{refreshMeta}")
		}
		$opts{refreshJaccard}=$opts{refresh} . "_jaccard.txt";
		if(! -e $opts{refreshJaccard}){
			usage("ERROR : jaccard data file for refresh not found, $opts{refreshJaccard}")
		}
		
	}
}

sub usage{
	
	print "\ngenerate_saucr_data.pl : generation of meatadata and jaccard scores for saucr\n";
	print "\nperl generate_saucr_data.pl [options]\n";
	print "Options are as follows:\n";
	print "\t--fpr.  The file provenence report to use, defaults to the system current fpr\n";
	print "\t--prefix. The prefix to use on output files.  Defaults to 'saucr'.  files are prefix_metadata.txt and prefix_jaccard.txt\n";
	print "\t--refresh. The prefix to use on previous files, for refresh.  If indicated, jaccard scores will only be calculated for new fingerprints\n";
	print "\t--help displays this usage message.\n";
		

	die "\n@_\n\n";
}

