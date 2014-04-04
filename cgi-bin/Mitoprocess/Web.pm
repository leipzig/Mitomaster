package Mitoprocess::Web;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::Parsing;
use Mitoprocess::Results;
use Mitoprocess::MitoResults;
use Mitoprocess::Genbank;
#use strict;
#use warnings;

use Bio::Tools::Run::StandAloneBlast;
use IO::String;
use File::Basename;
use List::Util;
use CGI qw(:standard :session);
use HTML::Table::FromDatabase;
use Carp;

#for formatting gb and hap fracs
use Number::Format qw(:subs :vars);
$DECIMAL_DIGITS = 2;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT =
  qw(print_conservation_popup process_form process_post print_checkbox_popup display_sequence_form display_genbank_form display_snvs_form print_summary_preload_html print_summary_html display_web_page footer print_details_html print_polymorphisms_coding_html);


sub process_genbank_ids {
    my $query   = shift;
    my $genbank = $query->param('seqInputName');

    $genbank =~ s/\s//g;

    #a range
    my @terms;
    my $term;

    if ( $genbank =~ m/([A-Za-z]+\d+\.?\d?)\-([A-Za-z]+\d+\.?\d?)/ ) {

        my $fromGb = $1;
        my $toGb   = $2;

        #get numeric component
        #they are using 5.8.8 so no named captures
        $fromGb =~ m/([A-Za-z]+)(\d+)(\.\d)?/;
        my $sprefix  = $1;
        my $startNum = $2;
        $toGb =~ m/([A-Za-z]+)(\d+)(\.\d)?/;
        my $eprefix = $1;
        my $endNum  = $2;
        if ( $sprefix ne $eprefix ) {
            $query->append(
                -name => 'error_message',
                -values =>
                  ["Start ($sprefix) and end prefix ($eprefix) must match."]
            );
        }
        for ( my $gbi = $startNum ; $gbi <= $endNum ; $gbi++ ) {
            push @terms, $sprefix . $gbi;
        }
        $term = join ',', @terms;
    }
    else {
        $term = $query->param('seqInputName');
    }
    my $gb = Mitoprocess::Genbank->new( query => $term );
    my $sequences = $gb->get_fastas();
    return \$sequences;
}

sub process_post {
    my $query = shift;
    my $run_id;
    if ( validate_form($query) ) {

        # checking the attached file

        # file size limit
        $CGI::POST_MAX = 1024 * 5000;

        #getting fasta file path
        #old name

        my $filename = $query->param("file");

 #defining the sequence hash to be used later (this hash stores input sequences)
        my $seqfasta = _get_entries( "file", $filename, $query );
        return $seqfasta;
    }
    else {
        print "invalid query\n";
    }
}


sub process_form {
    my $query = shift;
    my $run_id;
    if ( validate_form($query) ) {

        # checking the attached file

        # file size limit
        $CGI::POST_MAX = 1024 * 5000;

        my $containerType = "";

        #sequence count
        my $totalseqcount = 0;

        #reading the sequence from the posted webform
        my $textAreaContents = $query->param('seqInputName');

        my $filename = $query->param("file");

        #setting sequence type
        if ($filename) {
            $containerType = "file";
        }
        else {
            $containerType = "single";
        }

        print "selected type: $containerType filename: $filename<br/>\n"
          if (DEBUG);

 #defining the sequence hash to be used later (this hash stores input sequences)
        my $seqfasta = _get_entries( $containerType, $filename, $query );
        return $seqfasta;
    }
    return;
}


sub validate_form {
    my $query = shift;

    if ( !$query->param("file") or ($query->param("file") eq "") ) {
        if ( $query->param('fileType') eq 'sequences'
            && length( $query->param('seqInputName') ) < MIN_SEQ_LENGTH )
        {
            $query->append(
                -name   => 'error_message',
                -values => ["Sequence length less than "
                  . MIN_SEQ_LENGTH
                  . "<br/>\n"]
            );

            #display_sequence_form( $query );
            return 0;
        }
        elsif ( $query->param('fileType') eq 'genbank' ) {

            #http://www.ncbi.nlm.nih.gov/Sequin/acc.html
            #1 letter + 5 numerals OR 2 letters + 6 numerals
            #GI Number (GenInfo Identifier): any series of digits
            unless ( $query->param('seqInputName') =~
                /^([A-Za-z]{1,2}\d{5,6}\.?\d?)\-([A-Za-z]{1,2}\d{5,6}\.?\d?)$/
                || $query->param('seqInputName') =~
                /^(([A-Za-z]{1,2}\d{5,6}\.?\d?|[0-9]+),?\s?)+$/ )
            {
                $query->append(
                    -name => 'error_message',
                    -values =>
["The GenBank identifier(s) do not appear to be an accession, accession.version or GI ids"
                      . "<br/>\n"]
                );
                return 0;
            }
        }
    }

    # Form OK - return success
    return 1;
}



sub get_head {
    my $page_title = shift;
    my $WEB_HOME  = WEB_HOME;
    my $SITE_NAME = SITE_NAME . " " . $page_title;
    print <<END_HTML;
    <head>
        <meta charset="utf-8">
        <title>$SITE_NAME</title>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta name="description" content="">
        <meta name="author" content="">

        <!-- Le styles -->
        <link href="$WEB_HOME/css/styles.css" rel="stylesheet">
        <link href="$WEB_HOME/bootstrap/css/bootstrap.css" rel="stylesheet">
        <link href="$WEB_HOME/bootstrap/css/bootstrap-responsive.css" rel="stylesheet">
        

        
        <!--Theme-->
        <link href="$WEB_HOME/jquery-ui-bootstrap/css/custom-theme/jquery-ui-1.8.16.custom.css" rel="stylesheet">

        <!--datatables-->
        <link href="$WEB_HOME/css/dataTables.bootstrap.css" rel="stylesheet">
        <link href="$WEB_HOME/TableTools/media/css/TableTools.css" rel="stylesheet">

END_HTML
}

sub display_web_page {
    my ($refp) = shift;
    print <<END_HTML;
        <!DOCTYPE html>
        <html lang="en">
END_HTML
    my %pages = (
        'sequences'  => 'Sequence Query',
        'genbank'    => 'GenBank Query',
        'snvs'       => 'SNV Query',
        'species'       => 'Conservation Species',
        'webservice' => 'Webservices'
    );
    get_head( $pages{$refp} );
    print "</head>";
    nav_bar($refp);
}

sub nav_bar {
    my ($refp) = shift;
    my $SITE_NAME = SITE_NAME;
    my %pages =
      ( 'sequences' => '', 'genbank' => '', 'snvs' => '', 'species' => '', 'webservice' => '' );
    $pages{$refp} = 'class="active"';
    print <<END_HTML;
    <body data-spy="scroll" data-target=".subnav" data-offset="50">
    <div class="navbar navbar-fixed-top">
        <div class="navbar-inner">
          <div class="container">
            <a class="btn btn-navbar" data-toggle="collapse" data-target=".nav-collapse">
            <span class="icon-bar"></span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
            </a>
            <a class="brand" href="index.cgi">$SITE_NAME</a>
            <div class="nav-collapse">
              <ul class="nav">
                <li $pages{sequences}><a href="index.cgi">Sequence Query</a></li>
                <li $pages{genbank}><a href="index_genbank.cgi">GenBank Query</a></li>
                <li $pages{snvs}><a href="index_snvs.cgi">SNV Query</a></li>
                <li $pages{species}><a href="cons_checkboxes.cgi">Species</a></li>
                <li $pages{webservice}><a href="webservice.cgi">Use Web API</a></li>
              </ul>
              <ul class="nav pull-right">
                <li><a href="http://mitomaster.research.chop.edu/bin/view.pl/MITOMASTER/MitomasterInstructions">Help</a></li>
                <li class="divider-vertical"></li>
        
              </ul>
            </div><!-- /.nav-collapse -->
          </div>
        </div><!-- /navbar-inner -->
      </div>
<!--added for footer benefit -->
		  <div class="content">
	<div class="wrapper">
	 <div class="proper-content">
END_HTML
}

sub footer {
	print <<END_HTML;
	</div><!-- /.proper-content -->
    <div class="push"></div>
    </div><!-- /.wrapper -->
    <div class="footer-wrapper">
        <footer>
            <div id="chopbar-container">
                  <div id="chopbar">
                      <p>
                        <a class="cbmi" href="http://www.research.chop.edu/programs/cbmi/">Center for Biomedical Informatics</a> 
                        <a class="chop-research" href="http://www.research.chop.edu/">The Children's Hospital of Philadelphia | Research Institute</a>
                      </p>   
                  </div>
            </div>
        </footer>
    </div>
    </div>
END_HTML
}


sub display_sequence_form {
    my $query = shift;

    # Remove any potentially malicious HTML tags
    if($query->param('fasta_path')){
        my $tmp_fp =~ $query->param('fasta_path') || undef;
        $tmp_fp =~ s/<([^>]|\n)*>//g;
        $query->param( -name => 'fasta_path', -values => $tmp_fp );
    }
    #print "$fasta_path - $sequence_text - $containerType <br/>\n" if(DEBUG);
    # Build "selected" HTML for the "You are submitting" radio buttons
    my $WEB_HOME = WEB_HOME;

    display_web_page('sequences');

    # Display the form
    print <<END_HTML;
<div class="container">
<form id="seqAnalysisForm" name="seqAnalysisForm" action="index.cgi" method="POST" enctype="multipart/form-data">
<input type="hidden" name="fileType" value="sequences">
   <h1>Submit sequence <small>as file upload or paste into text area</small></h1>
</div>
	<div id="consRow">
		<div class="container">
			<div class="row">
					<div class="span2">
					<h2>Step 1</h2>
					</div>
					<div class="span10">
					<h3><a href="cons_checkboxes.cgi">Select conservation species of interest</a></h3>
					</div>
			</div>
		</div>
	</div>
	<div class="choiceRow">
	<div class="container">
<div class="row">
<div class="span2">
<h2>Step 2<br/>&nbsp;Option 1</h2>
</div>
<div class="span10">    



END_HTML
    if ( $query->param(-name => 'error_message') ) {
        my @errors = $query->param( -name => 'error_message');
        print "<div class=\"alert alert-error\">
                <button type=\"button\" class=\"close\" data-dismiss=\"alert\">&times;</button>
                <strong>Error:</strong> "
          . join('<br/>',@errors) . "</div>";
    }
    print <<END_HTML;
<h3>Fasta file <small>use this for more than one sequence</small></h3>



<div class="fileupload fileupload-new" data-provides="fileupload">
  <span class="btn btn-file"><span class="fileupload-new">Select file</span><span class="fileupload-exists">Change</span><input type="file" name="file" id="file"/></span>
  <span class="fileupload-preview"></span>
  <a href="#" class="close fileupload-exists" data-dismiss="fileupload" style="float: none">&#215;</a>
</div>

<button type="submit" value="Submit" name="submit" class="btn btn-primary">Submit</button>
</div>
</div>
</div>
</div>
<div class="choiceRow">
<div class="container">
<div class="row">
<div class="span2">
<h2>&nbsp;Option 2</h2>
</div>
<div class="span10">
<h3>Copy-paste</h3>
<h4><small></small></h4>
END_HTML
    
    print $query->textarea(
        -name    => 'seqInputName',
        -class      => 'sequence-textarea',
        -default => '',
        -rows    => 10,
        -columns => 100,
		-placeholder => "Paste a single raw unheadered sequence here. Minimum length:" . MIN_SEQ_LENGTH . "bp"
    );
    print <<END_HTML;
<p><button type="submit" value="Submit" name="submit" class="btn btn-primary">Submit</button></p>
</div><!--span10--> 
</div><!-- /.row -->
</form>
</div><!--container-->
</div><!--choicerow-->
END_HTML
footer();
    _javascript();

    _statcounter(); #this is the only page we track
    print <<END_HTML;
</body>
</html>
END_HTML

}

sub display_genbank_form {
    my $query = shift;

    # Remove any potentially malicious HTML tags
    if($query->param('fasta_path')){
        my $tmp_fp =~ $query->param('fasta_path');
        $tmp_fp =~ s/<([^>]|\n)*>//g;
        $query->param( -name => 'fasta_path', -values => $tmp_fp );
    }
    #print "$fasta_path - $sequence_text - $containerType <br/>\n" if(DEBUG);
    # Build "selected" HTML for the "You are submitting" radio buttons
    my $WEB_HOME = WEB_HOME;

    display_web_page('genbank');

    # Display the form
    print <<END_HTML;


<div class="container">
   <h1>Submit GenBank Identifiers</h1>
</div>
	<div id="consRow">
		<div class="container">
			<div class="row">
					<div class="span2">
					<h2>Step 1</h2>
					</div>
					<div class="span10">
					<h3><a href="cons_checkboxes.cgi">Select conservation species of interest</a></h3>
					</div>
			</div>
		</div>
	</div>
	<div class="choiceRow">
	<div class="container">
    <form id="seqAnalysisForm" name="seqAnalysisForm" action="index_genbank.cgi" method="POST" enctype="multipart/form-data">
    <input type="hidden" name="fileType" value="genbank">
END_HTML
    if ( $query->param(-name => 'error_message') ) {
        print "<div class=\"alert alert-error\">
                <button type=\"button\" class=\"close\" data-dismiss=\"alert\">&times;</button>
                <strong>Error:</strong> "
          . $query->param(-name => 'error_message') . "</div>";
    }
    print <<END_HTML;
    <div class="row">
	<div class="span2">
	<h2>Step 2</h2>
	</div>
<div class="span10">
<h3>Input Accession(s) or GI(s)</h3>
<p>Allowed values: singletons (<a href="#" class="genbank_link">EF060316</a> <a href="#" class="genbank_link">93116889</a>), a comma-separated list (<a href="#" class="genbank_link">EF060316,302376313,DQ112752</a>), or a range (<a href="#" class="genbank_link">AF346963-AF346968</a>)</p>

<div class="input-append">
    <input type="text" name="seqInputName" id="gb_accession" placeholder="GenBank identifier"  class="span3 search-query">
    <button type="submit" value="Submit" name="submit" class="btn">Search</button>
  </div>
  
<!--see mitomap.js for manual activation-->
</div> </div><!-- /.row -->
</form>
</div>
</div><!--container-->


END_HTML
footer();
    _javascript();
    print <<END_HTML;
    <script src="$WEB_HOME/js/mitoids.js"></script>
    <script src="$WEB_HOME/js/genbank.js"></script>
</body>
</html>
END_HTML

}





sub display_snvs_form {
    my $query = shift;

    if($query->param('fasta_path')){
        my $tmp_fp =~ $query->param('fasta_path');
        $tmp_fp =~ s/<([^>]|\n)*>//g;
        $query->param( -name => 'fasta_path', -values => $tmp_fp );
    }
    #print "$fasta_path - $sequence_text - $containerType <br/>\n" if(DEBUG);
    # Build "selected" HTML for the "You are submitting" radio buttons
    my $WEB_HOME = WEB_HOME;

    display_web_page('snvs');

    print <<END_HTML;
    <div class="container">
       <h1>Submit variant list <small>as file upload or paste into text area</small></h1>

        <form id="seqAnalysisForm" name="seqAnalysisForm" action="index_snvs.cgi" method="POST" enctype="multipart/form-data">
        <input type="hidden" name="fileType" value="snvlist">
	</div>
	<div id="consRow">
		<div class="container">
			<div class="row">
					<div class="span2">
					<h2>Step 1</h2>
					</div>
					<div class="span10">
					<h3><a href="cons_checkboxes.cgi">Select conservation species of interest</a></h3>
					<label class="checkbox">
                      <input type="checkbox" name="haplotype" value="true">
                      Compute haplotype (only recommended if supply a full set of SNVs)
                    </label>
					</div>
			</div>
		</div>
	</div>
	
	<div class="choiceRow">
	<div class="container">
<div class="row">
<div class="span2">
<h2>Step 2<br/>&nbsp;Option 1</h2>
</div>
<div class="span4">
END_HTML
   if ( $query->param(-name=>'error_message') ) {
      my @errors = $query->param(-name=>'error_message');
        print "<div class=\"alert alert-error\">
            <button type=\"button\" class=\"close\" data-dismiss=\"alert\">&times;</button>
            <strong>Error:</strong> "
          . join('<br/>',@errors) ."</div>";
    }
    print <<END_HTML;


<h3>Copy-paste</h3>
END_HTML
    print $query->textarea(
        -name    => 'seqInputName',
        -id      => 'snv_list',
        -class      => 'snv-textarea',
        -default => '',
        -rows    => 10,
        -columns => 30,
		-placeholder => 'Paste your variants here - choose one of the formats to the right'
    );
    print <<END_HTML;
	<p><button type="submit" value="Submit" name="submit" class="btn btn-primary">Submit</button></p>

</div>
<div class="span6">
<h4>The variant list can follow either of these formats:</h4>
<ul>
    <li>Tabular (<a href="#" id="snv_tabular" class="snv_link" name="snv_tabular">load example</a>)
    </li>
    <li style="list-style: none">
        <ul>
            <li>Either no header or the following:
                <pre>
sample pos ref var
</pre>
            </li>
            <li>Four tab-delimited columns:
                <ul>
                    <li>sample_name (with no spaces)
                    </li>
                    <li>mitochondrial position
                    </li>
                    <li>reference base
                    </li>
                    <li>variant base
                    </li>
                </ul>
            </li>
            <li>Insertions are designated using a ':' or empty cell for the reference
            </li>
            <li>Deletions are designated using a ':' or empty cell for the variant
            </li>
            <li>Lines beginning in # are ignored for user comments
            </li>
        </ul>
    </li>
    <li>Compact (<a href="#" id="snv_compact" class="snv_link" name="snv_compact">load example</a>)
    </li>
    <li style="list-style: none">
        <ul>
            <li>One variant per line, may be tab-prefaced with a sample column
            </li>
            <li>Substitutions: 9028T or C9028T
            </li>
            <li>Insertions: C573CC or 573insC or 573.1C or :573C
            </li>
            <li>Deletions: 291d or A291d or A291del or A291- or A291:
            </li>
        </ul>
    </li>
</ul>

</div>
</div><!-- /.row -->
</div>
</div>
<div class="choiceRow">
    <div class="container">
        <div class="row">
            <div class="span2">
                <h2>
                    &nbsp;Option 2
                </h2>
            </div>
            <div class="span4">
                <h3>
                    Select a variant file<br>
                </h3>
                <div class="fileupload fileupload-new" data-provides="fileupload">
                    <span class="btn btn-file"><span class="fileupload-new">Select file</span><span class="fileupload-exists">Change</span><input type="file"></span>
                    <a href="#" class="close fileupload-exists" data-dismiss="fileupload" style="float: none">&#215;</a>
                </div>
                <p>
                    <button type="submit" value="Submit" name="submit" class="btn btn-primary">Submit</button>
                </p>
            </div>
            <div class="span6"></div>
        </div><!-- /.row -->
    </div><!-- choicerow-->
</div><!--container-->

</form>
END_HTML
footer();
    _javascript();
    print <<END_HTML;
</body>
</html>
END_HTML
}




sub species_checkboxes {
    my $query = shift;    
    my %checked_species;
    if($query->cookie("CGISESSID")){
        my $sid = $query->cookie("CGISESSID");
        my $session = new CGI::Session(undef, $sid, {Directory=>'/tmp'});
        if($session->param(-name=>'species')){
            my $specptr = $session->param(-name=>'species');
            foreach my $checked(@{$specptr}){
                $checked_species{$checked}=1;
            }
        }
    }
#my $mainheading="SELECT DISTINCT "
#    my $speciesSql =
#"SELECT  species,replace(species,'_',' ') as species_name, common_name from mitomasterb.species ORDER BY species_order";
    my @species;
    my %labels;
    my $dsn = getDSN();

    my $i;
    my @mammys = qw(Mammals Non-Mammals);
    foreach my $mammy (@mammys) {
        print "<h3>$mammy</h3>";
        my $mammalsSQL =
"SELECT subheading,replace(subheading,' ','_') as subheading_link, MIN(species_order) as suborder from mitomasterb.species WHERE mainheading=? GROUP BY subheading ORDER BY suborder, subheading;";
        my $sth = $dsn->prepare($mammalsSQL);
        $sth->execute($mammy);

        #does opening one collapse the others?
        my $accordion = '';
        $accordion = "data-parent=\"#accordion2\"" if (ACCORDION_TOGGLE);

        while ( my $hash_ref = $sth->fetchrow_hashref ) {
            my $speciesSql =
"SELECT  species,replace(species,'_',' ') as species_name, common_name from mitomasterb.species WHERE subheading=? ORDER BY species_order";
            my $speciesth = $dsn->prepare($speciesSql);
            $speciesth->execute( $hash_ref->{subheading} );
            print <<END_HTML;
<div class="accordion" id="accordion2">
    <div class="accordion-group">
      <div class="accordion-heading">
        <a class="accordion-toggle" data-toggle="collapse" $accordion href="#$hash_ref->{subheading_link}">
          $hash_ref->{subheading}
        </a>
        
      </div>

      <div id="$hash_ref->{subheading_link}" class="accordion-body collapse">
        <div class="accordion-inner">


<a class="btn btn-mini" rel="$hash_ref->{subheading_link}" href="#select_all">Select All</a>
<a class="btn btn-mini" rel="$hash_ref->{subheading_link}" href="#select_none">Select None</a>

<fieldset id="$hash_ref->{subheading_link}">
END_HTML
            print "<table class=\"table\">";
            $i = 0;
            while ( my $hash_ref = $speciesth->fetchrow_hashref ) {
                if ( $i == 0 ) { print "<tr>"; }
                print "<td width=\"25%\"><label class=\"checkbox inline\">";
                print
"<input type=\"checkbox\" class=\"species_check\" name=\"species\" id=\"$hash_ref->{'species'}\" value=\"$hash_ref->{'species'}\"";
                if(defined $checked_species{$hash_ref->{'species'}}){
                    print " checked=\"yes\"";
                }
                print ">";
                print
"<a href=\"#\" class=\"tooltip-species\" rel=\"tooltip\" title=\"$hash_ref->{'common_name'}\">";
                print $hash_ref->{'species_name'};
                print "</a></label></td>";
                $i++;
                if ( $i == 4 ) { print "</tr>"; $i = 0 }
            }

            unless ( $i == 0 ) {
                while ( $i < 4 ) { print "<td  width=\"25%\"></td>"; $i++ }
                print "</tr>";
            }
            print "</table>";

            print <<END_HTML;
</fieldset>      
                 </div><!--accordion-inner-->
      </div><!--accordion-body collapse-->
    </div><!--accordion-group-->
</div><!--accordion-->
END_HTML
        }
    }

}


sub print_checkbox_popup {
    my ( $query ) = @_;
    
    display_web_page('species');
    print <<END_HTML;
    <div class="container">
	<h3>Conservation index species<div id="success"></div></h3>
	<h4><small>All species are used by default. Select a subset of species for the conservation index calculation in genes, tRNAs, and rRNAs.</small></h4>
	<h4><small>Mouse over to view common name.</small></h4>
    <form id="species_checkbox">
    <div id="all_species">
    <a class="btn btn-mini"  href="#select_everything">Select All</a>
    <a class="btn btn-mini"  href="#select_crs">Clear All (except rCRS)</a>
    <a class="btn btn-mini"  href="#select_primates">Primates Only</a>
END_HTML
species_checkboxes($query);
    print <<END_HTML;
    </div>
    </form>
</div>
END_HTML
footer();
_javascript();
print <<END_HTML;
</body>
</html>
END_HTML
}
    
sub print_conservation_popup {
    my ( $rcrs, $locus, $pos , $res, $queryid) = @_;
    my $WEB_HOME = WEB_HOME;
    my $query_name = get_query_name($queryid);
    print <<END_HTML;
            <!DOCTYPE html>
            <html lang="en">
END_HTML
    get_head('popup');
    _javascript();
    print <<END_HTML;
         <script type="text/javascript" charset="utf-8">
 			\$(document).ready(function() {
                \$('#conspopup').html( '<table cellpadding="0" cellspacing="0" border="0" class="table table-condensed table-striped table-bordered" id="conspopupdt"></table>' );
 				var oTable = \$('#conspopupdt').dataTable( {
 				     "bPaginate": false,
                        "bLengthChange": false,
                        "bFilter": false,
                        "bSort": false,
                        "bInfo": false,
                        "bAutoWidth": false,
                        // l - Length changing
                        //                                                 f - Filtering input
                        //                                                 t - The table!
                        //                                                 i - Information
                        //                                                 p - Pagination
                        //                                                 r - pRocessing
                        //                                                 < and > - div elements
                        //                                                 <"class" and > - div with a class
                        //                                                 Examples: <"wrapper"flipt>, <lf<t>ip>
                        // T is tabletools
                        "sDom": "tr",
 					"bProcessing": true,
 					"sAjaxSource": "conservation.cgi?rcrs=$rcrs&locus=$locus&pos=$pos&res=$res",
 						"aoColumns": [
                            { "sTitle": "Species", "sWidth": "10px" },
                            {"sTitle": "-5", "sWidth":"10px", "sClass":"Center"},
                            {"sTitle": "-4", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "-3", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "-2", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "-1", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "0", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "+1", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "+2", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "+3", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "+4", "sWidth":"10px", "sClass":"center"},
                            {"sTitle": "+5", "sWidth":"10px", "sClass":"center"}
 						]
 				} );
 			} );
 		</script>
 		</head>
 		<body>
 		<div class="container">
 		<h4>$query_name conservation grid</h4>
 		<p>Conservation denotes the percentage of residues that match <strong>rCRS</strong> among selected species. This query sequence mutation is depicted in isolation - i.e. other mutations in your sequence are not shown - and flanked by rCRS reference for clarity.
        Loci are shown in their transcribed orientation.</p>
                <div id="conspopup">
                </div>
        </div>
        </body></html>
END_HTML
}

sub print_polymorphisms_coding_html {
    my $WEB_HOME = WEB_HOME;
    print <<END_HTML;
            <!DOCTYPE html>
            <html lang="en">
END_HTML
    get_head('popup');
    _javascript();
print <<END_HTML;
<script type="text/javascript" charset="utf-8">
\$(document).ready(function() {
    \$('#aDataTable').dataTable( {
        "bProcessing": true,
        "bServerSide": true,
        "sAjaxSource": "http://resmitod.research.chop.edu/cgi-bin/polymorphisms_coding_json.cgi",
        "bPaginate": true,
        "bLengthChange" : true,
        "bFilter": false,
        "bSort": true,
        "bInfo": true,
        "bAutoWidth": true,
        "iDisplayLength": 50,
        "iDisplayStart": 0,
END_HTML
#_print_datatable_details();
print <<END_HTML;
});
});
</script>
</head>
<body>
<div class="container">
    <table cellpadding="0" cellspacing="0" border="0" class="table table-condensed table-striped table-bordered" id="aDataTable">
      <thead>
        <tr>
          <th width="20%">Locus</th>
          <th width="20%">Position</th>
          <th width="25%">Change</th>
          <th width="25%">Codon Count</th>
          <th width="15%">Codon Position</th>
          <th width="15%">AA Change</th>
          <th width="15%">GB Count</th>
          <th width="15%">Refs</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td colspan="5" class="dataTables_empty">Loading data from server</td>
        </tr>
      </tbody>
      <tfoot>
        <tr>
          <th>Locus</th>
          <th>Position</th>
          <th>Change</th>
          <th>Codon Count</th>
          <th>Codon Position</th>
          <th>AA Change</th>
          <th>GB Count</th>
          <th>Refs</th>
        </tr>
      </tfoot>
    </table>
    </div>
    </body></html>
END_HTML
}

sub get_details_data {
    my $query_id  = shift;

    print <<END_HTML;
         <script type="text/javascript">
             var aDataSet = [
END_HTML
my $haplo_seq_cnt=_details_table_innards($query_id);
my $hapname=get_query_haplogroup($query_id);
my $hapcol = ($haplo_seq_cnt>0) ? "{ \"sTitle\" : \"Freq % in $hapname\" }," : "";
_details_table_footer($hapcol);
}

sub get_all_details_data {
    my $run_id  = shift;
    my @query_ids=Mitoprocess::Results::get_query_ids($run_id );
    my $MITOMAP_HOME_CONF = MITOMAP_HOME;
    my $WEB_HOME = WEB_HOME;
    print <<END_HTML;
         <script type="text/javascript">
             var aDataSet = [
END_HTML
    my $haplo_seq_cnt=0;
    foreach my $query_id(@query_ids){
        $haplo_seq_cnt+=_details_table_innards($query_id);
    }
    my $hapcol = ($haplo_seq_cnt>0) ? "{ \"sTitle\" : \"Freq % in haplo\" }," : "";
    _details_table_footer($hapcol);
}

sub _details_table_innards {
    my $queryID  = shift;
my $sth_ptr  = get_var_report($queryID);
my $genbank_seq_cnt=get_genbank_seq_cnt(); #all seqs in genbank
my $hapname=get_query_haplogroup($queryID);
my $haplo_seq_cnt=get_haplo_total($hapname); #all seqs in genbank with this haplogroup
#if there is no haplogroup then don't show

my $MITOMAP_HOME_CONF = MITOMAP_HOME;
my $WEB_HOME = WEB_HOME;
    while ( my $hash_ref = $$sth_ptr->fetchrow_hashref ) {
         my $refs = get_refs($hash_ref->{rtmutid},$hash_ref->{mmutid},$hash_ref->{polyid});
         my ($pos,$ref,$alt)=($hash_ref->{tpos},$hash_ref->{tnt},$hash_ref->{qnt});
         my $refstring = ($refs) ? "<a href='$MITOMAP_HOME_CONF/cgi-bin/print_ref_list?refs=$refs;title=$ref>$alt at $pos'>refs</a>" : "";
         
         my $hapfrac = ($haplo_seq_cnt>0) ? format_number($hash_ref->{hap_frac}/$haplo_seq_cnt*100,2,1) : "";
         
         my $hapstring = ($haplo_seq_cnt>0) ? "<a href='".$MITOMAP_HOME_CONF."/cgi-bin/index_mitomap.cgi?title=Polymorphism+$ref-$alt+at+rCRS+position+$pos&pos=$pos&ref=$ref&alt=$alt'>".$hapfrac."</a>\",\"" : "";
         
         my $gbstring = "<a href='".$MITOMAP_HOME_CONF."/cgi-bin/index_mitomap.cgi?title=Polymorphism+".$hash_ref->{tnt}."-".$hash_ref->{qnt}."+at+rCRS+position+".$hash_ref->{tpos}."&pos=$pos&ref=$ref&alt=$alt'>".$hash_ref->{'gb_frac'}." (".format_number($hash_ref->{'gb_frac'}/$genbank_seq_cnt*100,2,1)."%)</a>";
        
        print "[\""
          . $hash_ref->{query} . "\"," . "\""
          . $hash_ref->{tpos} . "\"," . "\""
          . $hash_ref->{qpos} . "\"," . "\""
          . $hash_ref->{tnt} . "\"," . "\""
          . $hash_ref->{qnt} . "\"," . "\"";
         print $hash_ref->{ntchange} . "\"," . "\""
          . $hash_ref->{calc_locus} . "\"," . "\""
          . $hash_ref->{calc_aachange} . "\",". "\""
           . $gbstring . "\"," . "\""
           . $hapstring 
            . $refstring."\","
          . "\"<a href='popup.cgi?";
         print "rcrs=".$hash_ref->{tpos}
          . "&locus=".$hash_ref->{locus_id} 
          . "&pos=".$hash_ref->{conservation_pos} 
          . "&res=".$hash_ref->{variant_res}
          . "&query=".$queryID
          . "'>";
          print $hash_ref->{conservation}
          . "</a>\"," . "\""
          . $hash_ref->{patientphenotype}
          . "\"],\n";
    }
    return $haplo_seq_cnt;
}

sub _details_table_footer {
my $hapcol=shift;
print <<END_HTML;
                           ];
\$(document).ready(function() {
                           				\$('#detailTable').html( '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered" id="aDataTable"></table>' );
END_HTML

#<table id="aDataTable"> lives inside <div  id="detailTable">
_print_datatable_settings();
print <<END_HTML;
                           				    "aoColumns": [
                           				                            { "sTitle": "Query" },
                                            						{ "sTitle": "rCRS Position" },
                                            						{ "sTitle": "Query Position" },
                                            						{ "sTitle": "rCRS NT" },
                                            						{ "sTitle": "Query NT" },
                                                                    { "sTitle": "Mut type" },
                                                                    { "sTitle": "Locus" },
                                                                    { "sTitle": "tx Effect" },
                                                                    { "sTitle" : "GB Freq" },
                                                                    $hapcol
                                                                    { "sTitle" : "Refs" },
                                                                    { "sTitle": "Conservation" },
                                                                    { "sTitle": "Patient Report"}
                                            					]
                                            				} );//dtTableJsVar
//dtTableJSVar.\$("a[rel=popover]").popover().click(function(e) {e.preventDefault();});
// http://stackoverflow.com/questions/8130069/load-bootstrap-js-popover-content-with-ajax
// table.$("a[rel=popover]").popover().click(function(e) {e.preventDefault();});});
  dtTableJSVar.\$('.withajaxpopover').bind('click',function(myclickevent){
          myclickevent.preventDefault();
          var el=\$(this);
          \$.get(el.attr('data-load'),function(d){
              el.popover({placement: 'left', content: d}).popover('show');
          });
  });
                                            			} );
    
                                                                                              			
</script>
END_HTML
}




sub get_table_data {
    my $refp     = shift;
    my $runID    = shift;
    my $sth      = get_summary($runID);
    my $WEB_HOME = WEB_HOME;
    
    my $query = new CGI;
    my $sid = $query->cookie("CGISESSID") || undef;
    my $session = new CGI::Session(undef, $sid, {Directory=>'/tmp'});
    my $print_haplotype=0;
    
    if($session->param(-name=>'haplotype')){
        my $happtr = $session->param(-name=>'haplotype');
        foreach my $checked(@{$happtr}){
            if($checked eq 'true'){$print_haplotype=1;}
        }
    }
    
    print <<END_HTML;
         <script type="text/javascript">
             var aDataSet = [
END_HTML
    my $output;
    while ( my $hash_ref = $$sth->fetchrow_hashref ) {
        #my $poslist = $hash_ref->{poslist};
        #$poslist =~ s/;$//;
        #$poslist =~ s/;/,/g;
        $output.=
"[\"<span style='white-space:nowrap;'><a href='view_details.cgi?refp=$refp&queryid="
          . $hash_ref->{id}
          . "&runid=$runID'><span style='float:left' class='ui-icon ui-icon-circle-zoomin'></span>"
          . $hash_ref->{label}
          . "</a>&nbsp;&nbsp;&nbsp;</span>\"," . "\"";
          if($print_haplotype){$output.=$hash_ref->{haplogroup} . "\"," . "\""}
          $output.=$hash_ref->{varcount} . "\"," . "\""
          # . $hash_ref->{polycount} . "\"," . "\""
          # . $hash_ref->{patvariants} . "\"," . "\""
          . get_variant_list($hash_ref->{id})
          . "\"],\n";
    }
    print substr($output,0,length($output)-2);
    #"sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
    print <<END_HTML;
                               ];
\$(document).ready(function() {
              				\$('#summaryTable').html( '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered" id="aDataTable"></table>' );
END_HTML
    _print_datatable_settings();
    # { "sTitle": "Population Variants"},
    #   { "sTitle": "Patient Variants"},
    print <<END_HTML;
                                "aoColumns": [
                                  { "sTitle": "Sequence"},
END_HTML
if($print_haplotype){print  "{ \"sTitle\": \"Predicted Haplogroup\" },";}
print <<END_HTML;
                                  { "sTitle": "Total Variants"},
                                   { "sTitle": "Variants", "sWidth":"50%"}
                               ]
                             });
     });
    </script>
END_HTML
}

sub print_summary_preload_html {
    my ( $refp, $run_id ) = @_;
    print <<END_HTML;
            <!DOCTYPE html>
            <html lang="en">
END_HTML

    get_head($run_id);
    _javascript();
    print "</head>";
    nav_bar($refp);    #also starts body
    print <<END_HTML;
                <div class="container">
                <!-- I create the div for the progress bar. -->
                <div id="progresstext"><h4>Processing sequences</h4></div>
                <div id="progressbar"></div>

                <!-- I initialize the progress bar. -->
                <script type="text/javascript">
                    var cpt = 0
                    \$("#progressbar").progressbar({
                        value: cpt
                    });
                </script>

                <!-- I flush the response so the browser show the progress bar before the page is completely loaded. -->


END_HTML
}



sub print_summary_html {
    my ( $refp, $runID ) = @_;
    print "<script type=\"text/javascript\">\n";
    print "\$(\"#progressbar\").progressbar(\"destroy\")\n";
    print "\$(\"#progresstext\").hide();\n";
    print "</script>\n";
    _change_page_title( SITE_NAME . "_" . $runID );
    my $gg_pic = _get_genomegraph($runID);

#my $table = HTML::Table::FromDatabase->new( -sth => $$sth , -id => 'demo-dom');
#$table->print;
    get_table_data( $refp, $runID );

    #wijdemo();

    print <<END_HTML;
        <div class="header">
        
            <h2>
                Alignment Summary</h2>
            <h3>rCRS track view</h3>
        </div>
        <img src="$gg_pic" class="img-polaroid">
        <p>
            <h3>Sequence alignment<small><br/>select a sequence to show details or <a href="view_details.cgi?refp=$refp&runid=$runID"> show details for all sequences</a></small></h3>
            
            <div id="summaryTable">
            </div>
        </p>
    </div>
END_HTML
footer();
    print <<END_HTML;
    </body>
    </html>
END_HTML
}

sub print_details_html {
    my ( $refp, $queryID, $runID ) = @_;
    print <<END_HTML;
            <!DOCTYPE html>
            <html lang="en">
END_HTML
    my $queryname;
    if(defined $runID && !defined $queryID){
        $queryname = "all sequences";
        get_head($queryname);
        _javascript();
        get_all_details_data($runID);
    }else{
        $queryname = get_query_name($queryID);
        get_head($queryname);
        _javascript();
        get_details_data($queryID);
    }
    print "</head>";

    nav_bar($refp);    #also starts body

    my %pages = (
        'sequences'  => 'index.cgi',
        'genbank'    => 'index_genbank.cgi',
        'snvs'       => 'index_snvs.cgi',
        'webservice' => 'webservice.cgi'
    );
    my %prettyName = (
        'sequences' => 'Sequence run',
        'genbank'   => 'Genbank run',
        'snvs'      => 'SNV Run'
    );
    my $summaryName = $prettyName{$refp};
    my $summaryPage = $pages{$refp};

    print <<END_HTML;

    <div class="container">
    <ul class="breadcrumb">
      <li><a href="$summaryPage?runid=$runID">$summaryName $runID</a> <span class="divider">/</span></li>
      <li class="active">$queryname</li>
    </ul>
        <div class="header">
            <h2>Alignment Details</h2>
        </div>
        <p>
            <div id="detailTable">
            </div>
            
        </p>
    </div>
END_HTML
    #_modal();
    print <<END_HTML;
    </body>
    </html>
END_HTML
}



sub _change_page_title {
    my $page_title = shift;
    print <<END_HTML;
    <script type="text/javascript">
        \$(document).ready(function() {

          document.title = '$page_title';

        });
      </script>
END_HTML
}
#used by process_post, process_form
sub _get_entries {
    my ( $containerType, $filename, $query ) = @_;
    my $entries;
    if ( $containerType eq "file" ) {    # fasta file selected

        if ( !$filename ) {
            exit;
        }

        # making the filename safe
        my ( $name, $path, $extension ) = fileparse( $filename, '\..*' );
        $filename = $name . $extension;

        #clean up the file name
        my $safe_filename_characters = "a-zA-Z0-9_.-";
        $filename =~ tr/ /_/;
        $filename =~ s/[^$safe_filename_characters]//g;
        if ( $filename =~ /^([$safe_filename_characters]+)$/ ) {
            $filename = $1;
        }
        else {
            die "Filename contains invalid characters";
        }

        #getting the file handle
        my $upload_filehandle = $query->upload("file");

        #saving the file
        open( UPLOADFILE, ">" . UPLOADDIR . "/$filename" ) or die "$!";
        binmode UPLOADFILE;
        while (<$upload_filehandle>) { print UPLOADFILE; }
        close UPLOADFILE;

        #read the fasta file in a hash
        if ( $query->param('fileType') eq 'sequences' ) {
            $entries = loadFasta( UPLOADDIR . "/$filename" );
        }
        elsif ( $query->param('fileType') eq 'snvlist' ) {
            $entries = loadSNVs( UPLOADDIR . "/$filename", $query);
        }
        else {
            carp( "unknown fileType" . $query->param('fileType') );
        }

        # fasta sequences are stored in entries hash

    }
    else    # single sequence from textbox
    {
        if ( $query->param('fileType') eq 'sequences' ) {

            #reading the sequence from the posted webform
            my $textAreaContents = $query->param('seqInputName');
            $textAreaContents =~ s/[^a-zA-Z]//g;
            $textAreaContents =  uc $textAreaContents;
            $textAreaContents =~ s/U/T/g;
            $textAreaContents =~ s/\n//g;
            $textAreaContents =~ s/\r//g;
            
            #making a fasta sequence out of the read sequence
            # my $subjectstring = ">subjectSeq\n" . $textAreaContents;
            my $seqlabel = 'User Sequence';
            if ( ( substr $textAreaContents, 0, 1 ) eq '>' ) {
                $seqlabel = substr $textAreaContents, 1,
                  index( $textAreaContents, "\n" ) - 1;
                $entries->{$seqlabel} = substr $textAreaContents,
                  index( $textAreaContents, "\n" ) + 1;
            }
            else {
                $entries->{$seqlabel} = $textAreaContents;
            }
        }
        elsif ( $query->param('fileType') eq 'genbank' ) {
            my $tmpFileName = UPLOADDIR . "/tmp_genbank_ids" . time;
            my $sequences   = process_genbank_ids($query);

            #saving the file
            open( UPLOADFILE, ">$tmpFileName" ) or die "$!";
            binmode UPLOADFILE;
            print UPLOADFILE ${$sequences};
            close UPLOADFILE;

            #read the fasta file in a hash
            $entries = loadFasta($tmpFileName);
        }
        elsif ( $query->param('fileType') eq 'snvlist' ) {

            #snvlist
            my $tmpFileName = UPLOADDIR . "/tmp_variant_list" . time;
            open( UPLOADFILE, ">$tmpFileName" ) or die "$!";

            #binmode UPLOADFILE;
            print UPLOADFILE $query->param('seqInputName') or die "$!";
            close UPLOADFILE;

            #read the snvs in a hash
            $entries = loadSNVs($tmpFileName,$query);
        }
        else {
            $query->append(
                -name   => 'error_message',
                -values => ["Unknown file type."]
            );
        }
    }    #textbox
    return $entries;
}
sub _get_genomegraph {
    my $run_id      = shift;
    my $time        = time;
    my $tmpFileName = GFFDIR . "/" . $run_id . "_" . $time;
    my $gff3        = rCRS_gff3 . get_summary_gff3($run_id);
    open( UPLOADFILE, ">" . $tmpFileName ) or die "$!";
    print UPLOADFILE $gff3 or die "$!";
    close UPLOADFILE;
    my $outfile = TEMP_PICS . "/$time.png";
    my $cmd =
        "cat $tmpFileName | "
      . GENOMETOOLS_HOME
      . "/bin/gt sketch "
      . TEMP_PICS_PATH
      . "$time.png";
    system($cmd);
    return $outfile;
}
sub _javascript {
    my $WEB_HOME = WEB_HOME;
    print <<END_HTML;
    <script type="text/javascript" src="http://platform.twitter.com/widgets.js"></script>
    <script src="$WEB_HOME/js/jquery-1.7.2.min.js"></script>
    <script src="$WEB_HOME/js/jquery-ui-1.8.20.custom.min.js"></script>
    <script src="$WEB_HOME/bootstrap/js/bootstrap.js"></script>
    <script src="$WEB_HOME/js/mitomap.js"></script>
     <script src="$WEB_HOME/js/jquery.dataTables.js"></script>
     <script src="$WEB_HOME/js/dataTables.bootstrap.js"></script>
    <script src="$WEB_HOME/js/ZeroClipboard.js"></script>
    <script src="$WEB_HOME/js/TableTools.js"></script>
      <script src="$WEB_HOME/js/spin.min.js"></script> 
END_HTML
}

sub _modal {
    print <<END_HTML;
    <!-- Modal -->
    <div id="myModal" class="modal hide  wide-modal">
    <div class="modal-header">
    Conservation denotes the percentage of residues that match <strong>rCRS</strong> among selected species. This query sequence mutation is depicted in isolation - i.e. other mutations in your sequence are not shown - and flanked by rCRS reference for clarity.
    Loci are shown in their transcribed orientation.
    </div>
      <div class="modal-body">
      </div>
      <div class="modal-footer">
        
        <button class="btn" data-backdrop="false" data-dismiss="modal" aria-hidden="true">Close</button>
      </div>
    </div>
END_HTML
}

sub _print_datatable_details {
    my $WEB_HOME = WEB_HOME;
        print <<END_HTML;
    // l - Length changing
    //                                                 f - Filtering input
    //                                                 t - The table!
    //                                                 i - Information
    //                                                 p - Pagination
    //                                                 r - pRocessing
    //                                                 < and > - div elements
    //                                                 <"class" and > - div with a class
    //                                                 Examples: <"wrapper"flipt>, <lf<t>ip>
    // T is tabletools
    "sDom": "<'row'<'span6'l><'span6'<'tabletoolsclass'T>>>t<'row'<'span6'i><'span6'p>>",
    "oTableTools": {
                "sSwfPath": "$WEB_HOME/TableTools/media/swf/copy_csv_xls_pdf.swf",
                "aButtons": [
                "copy",
                {
                    "sExtends": "xls",
                    "sFileName": "*.tab",
                    "sToolTip": "Save as tab-delimited"
                },
                    {
                        "sExtends": "csv",
                        "sFileName": "*.csv",
                    },
                    {
                        "sExtends": "pdf",
                        "sFileName": "*.pdf",
                    },
                    "print"
                ]
    		},
	"sPaginationType": "bootstrap",
	"oLanguage": {
             "sLengthMenu": 'Display <select>'+
                         '<option value="10">10</option>'+
                         '<option value="25">25</option>'+
                         '<option value="50">50</option>'+
                         '<option value="100">100</option>'+
                         '<option value="-1">All</option>'+
                         '</select> records'
           },
END_HTML
}

sub _print_datatable_settings {
    my $WEB_HOME = WEB_HOME;
    print <<END_HTML;
    var dtTableJSVar = \$('#aDataTable').dataTable( {
	    "bPaginate": true,
        "bLengthChange": true,
        "bFilter": false,
        "bSort": true,
        "bInfo": true,
        "bAutoWidth": true,
	    "aaData": aDataSet,
        
END_HTML
_print_datatable_details();
}


sub _statcounter {
print <<END_HTML;
    <!-- Start of StatCounter Code for Default Guide -->
    <script type="text/javascript">
    var sc_project=8389746; 
    var sc_invisible=1; 
    var sc_security="4b3e1143"; 
    </script>
    <script type="text/javascript"
    src="http://www.statcounter.com/counter/counter.js"></script>
    <noscript><div class="statcounter"><a title="hit counter"
    href="http://statcounter.com/free-hit-counter/"
    target="_blank"><img class="statcounter"
    src="http://c.statcounter.com/8389746/0/4b3e1143/1/"
    alt="hit counter"></a></div></noscript>
    <!-- End of StatCounter Code for Default Guide -->
END_HTML
}
1;
