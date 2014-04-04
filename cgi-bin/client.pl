#!/usr/bin/perl
use LWP::UserAgent;
use HTTP::Request::Common;

my $userAgent = LWP::UserAgent->new(timeout => 1800); #a half-hour

my $request = POST 'http://resmitod.research.chop.edu/mitomasterdev/websrvc.cgi',
    Content_Type => 'multipart/form-data',
    Content => [ file => ['test/mitotestnew.fasta'], fileType => 'sequences'];

my $response = $userAgent->request($request);
print $response->error_as_HTML . "\n" if $response->is_error;

if ($response->is_success) {
     print $response->decoded_content;
} else {
     die $response->status_line;
}
