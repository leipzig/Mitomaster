#!/usr/bin/perl
use LWP::UserAgent;
use HTTP::Request::Common;

use Benchmark;

my $start = new Benchmark;


my $userAgent = LWP::UserAgent->new(timeout => 1800); #a half-hour

my $request = POST 'http://mitomaster.mitomap.org/cgi-bin/websrvc.cgi',
    Content_Type => 'multipart/form-data',
    Content => [ file => [$ARGV[0]], fileType => 'snvlist'];

my $response = $userAgent->request($request);
print $response->error_as_HTML . "\n" if $response->is_error;

if ($response->is_success) {
     print $response->decoded_content;
} else {
     die $response->status_line;
}

my $end = new Benchmark;

my $elapsed = timediff ($end, $start);
print STDERR "Elapsed time: ", timestr ($elapsed), "\n";
