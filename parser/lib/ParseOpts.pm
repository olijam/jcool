package ParseOpts;
require Exporter;

@EXPORT_OK = qw/%Opts/;
@EXPORT_TAGS = (ALL => qw/%Opts/);
@ISA = qw/Exporter/;

BEGIN {
  for(my $i=0;$i<=$#::ARGV;$i++) {
    if($::ARGV[$i] =~ /^-{1,2}([^=]*)(=(.*))?/) {
      $Opts{$1} = $3 || 1;
      splice(@::ARGV,$i--,1);
    }
  }
}

1;
