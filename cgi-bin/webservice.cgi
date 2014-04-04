#!/usr/bin/perl
use CGI;
use CGI::Carp qw(fatalsToBrowser);

#use strict;
use warnings;
use Bio::Tools::Run::StandAloneBlast;
use IO::String;

use File::Basename;
use List::Util;

use Mitoprocess::Parsing;
use Mitoprocess::FormHandling;
use Mitoprocess::Results;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::Web;

# Create the CGI object
my $query = new CGI;

# Output the HTTP header
print $query->header();

# BLAST
$ENV{BLASTDIR}='/usr/local/bin/blast-2.2.25/bin/';

#######

# Process form if submitted; otherwise display it

display_web_page('webservice');
print <<END_HTML;
<div class="container">
<h1>Submitting through the web service</h1>
<p>Multiple sequences can be submitted in bulk using a simple POST mechanism. The results are a tab-delimited table, very similar to the detail page shown in the Mitomaster web application.</p>

The following fields are required:
<ul><li>file - a file upload</li><li>fileType - set this to "sequences" or "snvlist" depending on the format. See <a href="index_snvs.cgi">SNV Query tab</a> for accepted formats.</li></ul>
<h3>Perl client example</h3>
<pre>
#!/usr/bin/perl
use LWP::UserAgent;
use HTTP::Request::Common;

my \$userAgent = LWP::UserAgent->new(timeout => 1800); #a half-hour

#fileType can be sequences or snvlist
my \$request = POST 'http://mitomaster.mitomap.org/cgi-bin/websrvc.cgi',
    Content_Type => 'multipart/form-data',
    Content => [ file => ['mySequences.fasta'], fileType => 'sequences'];

my \$response = \$userAgent->request(\$request);
print \$response->error_as_HTML . "\n" if \$response->is_error;

if (\$response->is_success) {
     print \$response->decoded_content;
} else {
     die \$response->status_line;
}
</pre>


<h3>Python client example</h3>
<pre>
#!/usr/bin/python
from poster.encode import multipart_encode
from poster.streaminghttp import register_openers
import urllib2

register_openers()

#fileType can be sequences or snvlist
datagen, headers = multipart_encode({"file": open("mySequences.fasta"),'fileType':'sequences'})

request = urllib2.Request("http://mitomaster.mitomap.org/cgi-bin/websrvc.cgi", datagen, headers)
try:
     print urllib2.urlopen(request).read()
except urllib2.HTTPError, e:
     print "HTTP error: %d" % e.code
except urllib2.URLError, e:
     print "Network error: %s" % e.reason.args[1]
</pre>
</div><!--container-->
END_HTML
footer();
    print <<END_HTML;
</body>
</html>
END_HTML
