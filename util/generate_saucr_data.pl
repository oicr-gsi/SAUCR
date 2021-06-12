use strict;
use warnings;
use JSON::PP;
use Data::Dumper;
use File::Basename;
use Getopt::Long;


my %opts = ();
GetOptions(
	"fpr|f:s"        		=> \$opts{"fpr"},
	"refresh|r:s"			=> \$opts{"refresh"},
	"prefix|p:s"     		=> \$opts{"prefix"},
	"maxrecs|m:i"			=> \$opts{"maxrecs"},
	"help|h"			=> \$opts{"help"},
);
validate_options(\%opts);


(open my $META,">",$opts{prefix} . "_metadata.txt") || die "unable to open metadata file for output";
my @fields=qw/limskey project identity external_name library tissue_origin tissue_type group_id library_source_template_type prep_kit run platform pool file wfid wfrunid sites sites_nocovered sites_covered sites_wildtype sites_variant/;
print $META join("\t",@fields) . "\n";

(open my $JACCARD,">",$opts{prefix} ."_jaccard.txt") || die "unable to open jaccard file for output";


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




### collection of fin files for jaccard calcluations
my %fin;

### temp, for testing, restrict number of records
if($opts{maxrecs}){
	@recs=@recs[0..$opts{maxrecs}];
}

my %stats;

for my $rec(@recs){
	
	### store record information in a hash
	my @val=split /\t/,$rec;
	my %h;
	@h{@keys}=@val;
	
	### bypass if 1. no fin file or 2. no file
	next unless(-e $h{file});
	next unless($h{file}=~/fin$/);
	
	$stats{count}++;
	
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

print STDERR "starting jaccard score generation\n";
# module load sample-fingerprinting
#jaccard_coeff_pair

$stats{pairs}=($stats{count} * ($stats{count}+1))/2;

print STDERR "metadata has $stats{count} records.  Expecting $stats{pairs} single direction pairs\n";



### check records in the previous metadata report to asses which wfrundis are no longer used, and which are new
if($opts{refresh}){
	
	print STDERR "refreshing from $opts{refresh}\n";
	my @recs=`cat $opts{refreshMeta}`;chomp @recs;
	my @headers=split /\t/,shift @recs;
    #### update status in the fin hash based on the refresh.  
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

	(open my $LASTJACC,"<",$opts{refreshJaccard});
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


print STDERR "adding in new records\n";
### get any new fin files
my @new;
map{ push(@new,$_) if($fin{$_}{status} eq "new") }keys %fin;


### for all new run ids, generated jaccard scores againts ALL other ids, but only in one direction, alphabetically, by limskey
for my $runid1(@new){
	my $lk1=$fin{$runid1}{limskey};
	my $status1=$fin{$runid1}{status};
	my $fin1=$fin{$runid1}{fin};
	for my $runid2(sort keys %fin){
		my $lk2=$fin{$runid2}{limskey};
		my $status2=$fin{$runid2}{status};
		my $fin2=$fin{$runid2}{fin};
		
		my $check=0;
		my @lk_array=sort($lk1,$lk2);
		my @wf_array=sort($runid1,$runid2);
		
		if($status2 eq "new"){
			### if both are new, then we only want to generate the record if they are in alphabetical order, by limskey.  This will prevent generating the record twice
			$check=1 if($lk_array[0] eq $lk1);
		}else{
			$check=1;
			## if the second record is not new, generated the jaccard metric, and save.  the other direction won't be run
		}
		
		if($check){
			my $json=`perl $opts{script} -id1 $lk1 -fin1 $fin1 -id2 $lk2 -fin2 $fin2`;
			my $h=decode_json $json;
			my $jaccard=$$h{jaccard_score} || 0;
			my $covered=$$h{pair}{covered} || 0;
			my $lk_pair=join("\t",@lk_array);
			my $wfid_pair=join("_",@wf_array);
			
			print $JACCARD "$lk_pair\t$jaccard\t$covered\t$wfid_pair\n";
		}
		
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
	
	if($opts{refresh}){
		$opts{refreshMeta}=$opts{refresh} . "_metadata.txt";
		if(! -e $opts{refreshMeta}){
			usage("ERROR : meta data file for refresh not found, $opts{refreshMeta}")
		}
		$opts{refreshJaccard}=$opts{refresh} . "_jaccard.txt";
		if(! -e $opts{refreshJaccard}){
			usage("ERROR : jaccard data file for refresh not found, $opts{refreshJaccard}")
		}
		
	}
	
	if($opts{maxrecs}){
		print STDERR "limiting to $opts{maxrecs} fingerprintCollector records\n";
	}
	
	$opts{script}="/.mounts/labs/gsiprojects/gsi/jaccard/scripts/jaccard_coeff.pair.pl";
	
	
}

sub usage{
	
	print "\ngenerate_saucr_data.pl : generation of meatadata and jaccard scores for saucr\n";
	print "\nperl generate_saucr_data.pl [options]\n";
	print "Options are as follows:\n";
	print "\t--fpr.  The file provenence report to use, defaults to the system current fpr\n";
	print "\t--prefix. The prefix to use on output files.  Defaults to 'saucr'.  files are prefix_metadata.txt and prefix_jaccard.txt\n";
	print "\t--refresh. The prefix to use on previous files, for refresh.  If indicated, jaccard scores will only be calculated for new fingerprints\n";
	print "\t--maxrecs. The maxiumum number of fingerprint records to process from FPR, for testing\n";
	print "\t--help displays this usage message.\n";
		

	die "\n@_\n\n";
}

