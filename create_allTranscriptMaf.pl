#!/gsc/bin/perl

use strict;
use warnings;
use IO::File;
use Genome;

if( scalar( @ARGV ) != 2 ) {
  print STDERR "Usage: perl $0 <anno_file_dir> <output_maf>\n";
  exit 1;
}

# Location of annotation files to use, and the final merged maf file to create
my $anno_files_dir = $ARGV[0];
my $final_maf_file = $ARGV[1];

my %anno_files = map {chomp; m/files\/(.*).anno/; ($1, $_)} `ls $anno_files_dir/*.anno`;
my $out_maf_fh = IO::File->new( $final_maf_file, ">" ) or die "Cannot open $final_maf_file. $!";
$out_maf_fh->print( "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\t",
                    "End_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\t",
                    "Tumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tdbSNP_Val_Status\t",
                    "Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\t",
                    "Match_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\t",
                    "Match_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\t",
                    "Verification_Status\tValidation_Status\tMutation_Status\tSequencing_Phase\t",
                    "Sequence_Source\tValidation_Method\tScore\tBAM_File\tSequencer\t",
                    "chromosome_name\tstart\tstop\treference\tvariant\ttype\tgene_name\t",
                    "transcript_name\ttranscript_species\ttranscript_source\ttranscript_version\t",
                    "strand\ttranscript_status\ttrv_type\tc_position\tamino_acid_change\tucsc_cons\t",
                    "domain\tall_domains\tdeletion_substructures\ttranscript_error\n" );
foreach my $case ( keys %anno_files ) {
    my $anno_file = $anno_files{$case};
    my %enst=map{chomp; my @a=split(/\t/); ($_,1)}`cat $anno_file|cut -f8|sort -u|grep "^ENST"`;
    my $snv_anno_file = Genome::Sys->create_temp_file_path();
    my $indel_anno_file = Genome::Sys->create_temp_file_path();
    my $maf_file = Genome::Sys->create_temp_file_path();
    my $temp_maf_file = Genome::Sys->create_temp_file_path();

    # Separate the indels and snvs so that we can use the "gmt capture create-maf-file" tool
    my @anno_lines = `cat $anno_file`;
	foreach my $trid (keys %enst) {
		my $snv_anno_fh = IO::File->new( $snv_anno_file, ">" );
		my $indel_anno_fh = IO::File->new( $indel_anno_file, ">" );
		foreach my $line ( @anno_lines ) {
			my @anno_cols = split( /\t/, $line );
			#warn "$trid\n";
			$snv_anno_fh->print( $line ) if( $anno_cols[5] =~ m/SNP|DNP|TNP|ONP/ && $anno_cols[7] eq $trid);
			$indel_anno_fh->print( $line ) if( $anno_cols[5] =~ m/INS|DEL/ && $anno_cols[7] eq $trid);
		}
		$snv_anno_fh->close;
		$indel_anno_fh->close;
		#my @s=`cat $snv_anno_fh`;print "@s\n";
		my $maf_cmd = Genome::Model::Tools::Capture::CreateMafFile->create(
			snv_annotation_file => $snv_anno_file, snv_file => $snv_anno_file, output_file => $maf_file,
			indel_annotation_file => $indel_anno_file, indel_file => $indel_anno_file,
			somatic_status => "Somatic", genome_build => "37", platform => "-", center => "-",
			tumor_sample => $case, normal_sample => $case );
		( $maf_cmd->execute ) or die "Failed to run 'gmt capture create-maf-file' for variants in $case\n";
		$maf_cmd->delete;

		my @maf_lines = `cat $maf_file`;
		chomp( @maf_lines );
		foreach my $line ( @maf_lines ) {
			$out_maf_fh->print( "$line\n" ) unless( $line =~ m/^Hugo_Symbol/ );
			#print( "$line\n" ) unless( $line =~ m/^Hugo_Symbol/ );
		}
	}
}
$out_maf_fh->close;
