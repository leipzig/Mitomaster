package Mitoprocess::MitoWeb;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::Parsing;
use Mitoprocess::MitoResults;
use Mitoprocess::Web;

use strict;
use warnings;

use Bio::Tools::Run::StandAloneBlast;
use IO::String;
use File::Basename;
use List::Util;
use CGI qw(:standard);
use HTML::Table::FromDatabase;
use Carp;
use Number::Format qw(:subs :vars);
$DECIMAL_DIGITS = 2;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT =
  qw(print_summary_preload_html_mito print_summary_html_mito display_web_page_mito print_details_html_mito print_summary_html_haplogroup);
  
sub get_head {
    #<link href="/$WEB_HOME/bootstrap/js/google-code-prettify/prettify.css" rel="stylesheet">
    #<link href="http://cdn.wijmo.com/themes/rocket/jquery-wijmo.css" rel="stylesheet" type="text/css" title="rocket-jqueryui" />
    
    #<link href="/$WEB_HOME/bootstrap/css/bootstrap.css" rel="stylesheet">
    #<link href="/$WEB_HOME/bootstrap/css/bootstrap-responsive.css" rel="stylesheet">
    
    my $WEB_HOME=WEB_HOME;
    my $SITE_NAME=shift;
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

sub javascript {
    #<script src="/$WEB_HOME/bootstrap/js/google-code-prettify/prettify.js"></script>
    #<script src="/$WEB_HOME/jquery-ui-bootstrap/bootstrap/js/bootstrap-twipsy.js" rel="stylesheet">
    #<script src="/$WEB_HOME/jquery-ui-bootstrap/bootstrap/js/bootstrap-popover.js" rel="stylesheet">
    my $WEB_HOME=WEB_HOME;
    print <<END_HTML;
    <!-- Le javascript
    ================================================== -->
    <script type="text/javascript" src="http://platform.twitter.com/widgets.js"></script>
    <script src="$WEB_HOME/js/jquery-1.7.2.min.js"></script>
    <script src="$WEB_HOME/js/jquery-ui-1.8.20.custom.min.js"></script>
    <script src="$WEB_HOME/bootstrap/js/bootstrap.js"></script>
    <script src="$WEB_HOME/bootstrap/js/bootstrap.min.js"></script>
    <script src="$WEB_HOME/bootstrap/js/bootstrap.min.js"></script>
    <script src="$WEB_HOME/js/mitomap.js"></script>
     <script src="$WEB_HOME/js/jquery.dataTables.js"></script>
     <script src="$WEB_HOME/js/dataTables.bootstrap.js"></script>
    <script src="$WEB_HOME/js/ZeroClipboard.js"></script>
    <script src="$WEB_HOME/js/TableTools.js"></script>
    <!-- this is some old version of jquery that works with Mike's autocomplete function http://www.nodstrum.com/2007/09/19/autocompleter/
    <script type="text/javascript" src="$WEB_HOME/js/jquery-1.2.1.pack.js"></script> -->


END_HTML
}

sub nav_bar_s{
    print <<END_HTML;  
    <body data-spy="scroll" data-target=".subnav" data-offset="50">
END_HTML
}

sub get_genbank_data {
    my $queryID = shift;
    my $sth_ptr=get_genbank_report($queryID);
    my $WEB_HOME=WEB_HOME;
    my $sth1;
	my $sth2;
	my $sth3;
	my $haplogroup_total;
	my $haplogroup_cnt;
	my $genbank_cnt;
	my $genbank_frq;
	my $haplogroup_freq;
	my $genbank_seq_cnt=get_genbank_seq_cnt(); #all seqs in genbank
    print <<END_HTML;
         <script type="text/javascript">
             var aDataSet = [
END_HTML
 while ( my $hash_ref = $$sth_ptr->fetchrow_hashref ) {
     $hash_ref->{patientphenotype} =~ s/"//g; # dam'it. It contains " somewhere;
	 
	 $sth1=get_genbank_mut_haplofrq($hash_ref->{tpos},$hash_ref->{tnt},$hash_ref->{qnt},$hash_ref->{haplogroup});
     $haplogroup_cnt=$$sth1->fetchrow_hashref;
    
	 $sth2=get_genbank_mut_allfrq($hash_ref->{tpos},$hash_ref->{tnt},$hash_ref->{qnt});
	 
     $genbank_cnt=$$sth2->fetchrow_hashref;
	 $genbank_frq=format_number($genbank_cnt->{cnt}/$genbank_seq_cnt*100,2,1);
	 
	    $haplogroup_total=get_haplo_total($hash_ref->{haplogroup});
        
	   $haplogroup_freq=format_number($haplogroup_cnt->{cnt}/$haplogroup_total*100,2,1);
	   
     print "[\"".$hash_ref->{tpos} . "\","
       . "\"".$hash_ref->{qpos} . "\","
       . "\"".$hash_ref->{tnt} . "\","
       . "\"".$hash_ref->{qnt} . "\","
        . "\"".$hash_ref->{ntchange} . "\","
        . "\"".$hash_ref->{calc_locus} . "\","
           . "\"".$hash_ref->{cal_aachange} . "\","
		   . "\"".$haplogroup_freq. "\","
		   . "\"".$genbank_frq. "\","
                       . "\"<a href='#' rel='tooltip' title='Percent of selected species' data-content='content'>".$hash_ref->{conservation} . "</a>\","
                        . "\"".$hash_ref->{patientphenotype} . "\"],\n"; 
 }
 #          . "\"".$hash_ref->{polymorphism} . "\","
 #                               				    "sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
print <<END_HTML;  
                               ];
\$(document).ready(function() {
                               				\$('#detailTable').html( '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered" id="example"></table>' );
                               				\$('#example').dataTable( {
                               				    "bPaginate": true,
                                                "bLengthChange": true,
                                                "bFilter": false,
                                                "bSort": true,
                                                "bInfo": true,
                                                "bAutoWidth": true,
                               				    "aaData": aDataSet,
                               				    "sDom": "<'row'<'span6'l><'span6'<'tabletoolsclass'T>>>t<'row'<'span6'i><'span6'p>>",
                               				    "oTableTools": {
                               				                "sSwfPath": "$WEB_HOME/swf/copy_csv_xls_pdf.swf"
                                                		},
                                        		"sPaginationType": "bootstrap",
                                        	"oLanguage": {
                 "sLengthMenu": 'Display <select>'+
                             '<option value="10">10</option>'+
                             '<option value="25">25</option>'+
                             '<option value="50">50</option>'+
                             '<option value="-1">All</option>'+
                             '</select> records'
               },
                               				    "aoColumns": [
                                                						{ "sTitle": "rCRS Position" },
                                                						{ "sTitle": "Query Position" },
                                                						{ "sTitle": "rCRS NT" },
                                                						{ "sTitle": "Query NT" },
                                                                        { "sTitle": "Type of NT Change" },
                                                                        { "sTitle": "Locus" },
                                                                        { "sTitle": "AA Change" },
																		{ "sTitle": "Freq (%) in Haplogroup" },
																		{ "sTitle": "Freq (%) in GenBank" },
                                                                        { "sTitle": "Conservation" },
                                                                        { "sTitle": "Patient Report"}
                                                					]
                                                				} );	
                                                			} );
    </script>
END_HTML
#{ headerText: "Reported as Population Variant"},
}

sub get_table_data_genbank {
    my ($pos,$ref,$alt) = @_;
    my $sth=get_genbank($pos,$ref,$alt);
	my $sth1;
	my $sth2;
	my $haplogroup;
	my $haplogroup_total;
    my $WEB_HOME=WEB_HOME;
    print <<END_HTML;
         <script type="text/javascript">
             var aDataSet = [
	 
END_HTML
while ( my $hash_ref = $$sth->fetchrow_hashref ) {
# I will need to add st here
	     $sth1=get_haplo_cnt($pos,$ref,$alt,$hash_ref->{haplogroup});
         $haplogroup=$$sth1->fetchrow_hashref;
    
        my $reference_link='';
	    if($hash_ref->{reference}){
	        $reference_link="\"<span style='white-space:nowrap;'><a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&amp;cmd=Retrieve&amp;dopt=Abstract&amp;list_uids=".$hash_ref->{reference}."'>".$hash_ref->{reference} . "</a>&nbsp;&nbsp;&nbsp;</span>\","
	    }else{
	        $reference_link="\"\",";
	    }
        $haplogroup_total=get_haplo_total($hash_ref->{haplogroup});

     print "[\"<span style='white-space:nowrap;'><a href='http://www.ncbi.nlm.nih.gov/nuccore/".$hash_ref->{genbank_id}."'>".$hash_ref->{genbank_id} . "</a>&nbsp;&nbsp;&nbsp;</span>\","
        .$reference_link
        ."\"<a href='view_haplogroup.cgi?haplogroup=".$hash_ref->{haplogroup}."&pos=$pos&ref=$ref&alt=$alt"."'>".$hash_ref->{haplogroup} . "</a>\","
		. "\"".$haplogroup->{cnt} . "\","
		 . "\"".$haplogroup_total. "\","
		."\"<a href='view_mitomaster.cgi?genbankid=".$hash_ref->{genbank_id}."&haplogroup=".$hash_ref->{haplogroup}."&pos=$pos&ref=$ref&alt=$alt"."'>".'view'."</a>\"],\n";
 }
 #"sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
print <<END_HTML;  
                               ];
\$(document).ready(function() {
              				\$('#summaryTable').html( '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered" id="example"></table>' );
              				\$('#example').dataTable( {
	    "bPaginate": true,
        "bLengthChange": true,
        "bFilter": false,
        "bSort": true,
        "bInfo": true,
        "bAutoWidth": true,
	    "aaData": aDataSet,
"sDom": "<'row'<'span6'l><'span6'<'tabletoolsclass'T>>>t<'row'<'span6'i><'span6'p>>",
	    "oTableTools": {
	                "sSwfPath": "$WEB_HOME/TableTools/media/swf/copy_csv_xls_pdf.swf"
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
                                "aoColumns": [
                                  { "sTitle": "Genbank ID"}, 
                                  { "sTitle": "Pubmed Reference" },
                                  { "sTitle": "Predicted Haplogroup" },
								  { "sTitle": "#Haplogroup with given polymorphism " },
								  { "sTitle": "Total # Haplogroup in GenBank" },
                                  { "sTitle": "MitoMaster Results"}
                               ]
                             });
     });
    </script>
END_HTML
}


sub get_table_data_genbank_haplogroup {
    my ($pos,$ref,$alt,$haplo) = @_;
    my $sth=get_genbank_haplo($haplo);
	 
    my $WEB_HOME=WEB_HOME;
    print <<END_HTML;
         <script type="text/javascript">
             var aDataSet = [
	 
END_HTML
while ( my $hash_ref = $$sth->fetchrow_hashref ) {

     print "[\"<span style='white-space:nowrap;'><a href='http://www.ncbi.nlm.nih.gov/nuccore/".$hash_ref->{genbank_id}."'>".$hash_ref->{genbank_id} . "</a>&nbsp;&nbsp;&nbsp;</span>\","
        ."\"<span style='white-space:nowrap;'><a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&amp;cmd=Retrieve&amp;dopt=Abstract&amp;list_uids=".$hash_ref->{reference}."'>".$hash_ref->{reference} . "</a>&nbsp;&nbsp;&nbsp;</span>\","
		."\"<a href='view_mitomaster.cgi?genbankid=".$hash_ref->{genbank_id}."&haplogroup=$haplo&pos=$pos&ref=$ref&alt=$alt"."'>".'view'."</a>\"],\n";
 }

print <<END_HTML;  
                               ];
\$(document).ready(function() {
              				\$('#summaryTable').html( '<table cellpadding="0" cellspacing="0" border="0" class="table table-striped table-bordered" id="example"></table>' );
              				\$('#example').dataTable( {
	    "bPaginate": true,
        "bLengthChange": true,
        "bFilter": false,
        "bSort": true,
        "bInfo": true,
        "bAutoWidth": true,
	    "aaData": aDataSet,
"sDom": "<'row'<'span6'l><'span6'<'tabletoolsclass'T>>>t<'row'<'span6'i><'span6'p>>",
	    "oTableTools": {
	                "sSwfPath": "$WEB_HOME/TableTools/media/swf/copy_csv_xls_pdf.swf"
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
                                "aoColumns": [
                                  { "sTitle": "Genbank ID"}, 
                                  { "sTitle": "Pubmed Reference" },
                                  { "sTitle": "MitoMaster Results"}
                               ]
                             });
     });
    </script>
END_HTML
}


sub print_summary_preload_html_mito{
    my ($refp,$pagetitle) = @_;
     print <<END_HTML;
            <!DOCTYPE html>
            <html lang="en">
END_HTML

            get_head($pagetitle);
            javascript();
            print "</head>";
            nav_bar_s(); #also starts body
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

sub print_summary_html_mito {
my ($pos,$ref,$alt,$pagetitle) = @_;    
print "<script type=\"text/javascript\">\n";
print "\$(\"#progressbar\").progressbar(\"destroy\")\n";
print "\$(\"#progresstext\").hide();\n";
print "</script>\n";

     
    #my $table = HTML::Table::FromDatabase->new( -sth => $$sth , -id => 'demo-dom');
    #$table->print;
    get_table_data_genbank($pos,$ref,$alt);
     #wijdemo();
            print <<END_HTML;

        <div class="main demo">
            <h3>GenBank Record for $pagetitle <small><br/> click hyperlink to view GenBank Record, PubMed Reference, and Mitomaster running results</small></h3>
            <div id="summaryTable">
            </div>
        </div>
   
END_HTML
my $gb_seq_cnt=get_genbank_seq_cnt();  
    print <<END_HTML;

<h4><small>Results out of $gb_seq_cnt sequences in Mitomap's GenBank Set</small></h4>
 </div>
    </body>
    </html>
END_HTML
}


sub print_summary_html_haplogroup {
my ($pos,$ref,$alt,$haplo) = @_;    
print "<script type=\"text/javascript\">\n";
print "\$(\"#progressbar\").progressbar(\"destroy\")\n";
print "\$(\"#progresstext\").hide();\n";
print "</script>\n";
     
      get_table_data_genbank_haplogroup($pos,$ref,$alt,$haplo);
            print <<END_HTML;

        <div class="main demo">
            <h3>GenBank Record for haplogroup $haplo <small><br> click hyperlink to view GenBank Record, PubMed Reference, and Mitomaster running results</small></h3>
            <div id="summaryTable">
            </div>
        </div>
					<script>
	\$(function() {
		function log( message ) {
			\$( "<div/>" ).text( message ).prependTo( "#log" );
			\$( "#log" ).scrollTop( 0 );
		}

		\$( "#k" ).autocomplete({
			source: "search_haplogroup.php",
			minLength: 1,
			select: function( event, ui ) {
				log( ui.item ?
					"Selected: " + ui.item.value + " aka " + ui.item.id :
					"Nothing selected, input was " + this.value );
			}
		});
	});
	</script>
		<h3>Check other Haplogroup</h3>
		Haplogroup Name:  <input type="text" id="k" name="haplogroup" />
        <input type="button" id="Go" value="Search" onclick="location.href = 'http://resmitod.research.chop.edu/cgi-bin/view_haplogroup.cgi?pos=$pos&ref=$ref&alt=$alt&haplogroup='+document.getElementById('k').value;"/>

    </div>
END_HTML
    
    print <<END_HTML;
    </body>
    </html>
END_HTML
}


sub print_details_html_mito {
    my ($genbankid,$pos,$ref,$alt,$haplo) = @_; 
        print <<END_HTML;
            <!DOCTYPE html>
            <html lang="en">
END_HTML
    get_head($genbankid);
javascript();
        
    
    #my $table = HTML::Table::FromDatabase->new( -sth => $$sth , -id => 'demo-dom');
    #$table->print;
    get_genbank_data($genbankid);
     #wijdemo();
print "</head>";
nav_bar_s(); #also starts body
my %pages=('genbankid'=>'index_mitomap.cgi','genbank'=>'index_genbank.cgi','snvs'=>'index_snvs.cgi','webservice'=>'webservice.cgi');
my $summaryPage=$pages{'genbankid'};

#not sure why this is here if haplo is passed
#my $sth_haplo=get_haplo_id($genbankid);
#my $haplogroup_page_id=$$sth_haplo->fetchrow_hashref;
#    my $haplogroup_id=$haplogroup_page_id->{haplogroup_id};

print <<END_HTML;
    <div class="container">
        <div class="header">
            <h2>MitoMaster Running Results for $genbankid, haplogroup $haplo</h2>
        </div>
        <div class="main demo">
            <div id="detailTable">
            </div>
        </div>
		<br>
		<h3>Check other GenBank MitoMaster Results: <small><br>Using Accession IDs such as AY882410.1, etc </small><h3>
			<script>
	\$(function() {
		function log( message ) {
			\$( "<div/>" ).text( message ).prependTo( "#log" );
			\$( "#log" ).scrollTop( 0 );
		}

		\$( "#inputString" ).autocomplete({
			source: "search_genbank.php",
			minLength: 1,
			select: function( event, ui ) {
				log( ui.item ?
					"Selected: " + ui.item.value + " aka " + ui.item.id :
					"Nothing selected, input was " + this.value );
			}
		});
	});
	</script>
		GenBank ID:  <input type="text" value="" id="inputString" name="genbank" onkeyup="lookup(this.value);" onblur="fill();" />
        <input type="button" id="Go" value="Search" onclick="location.href = 'view_mitomaster.cgi?haplogroup=$haplo&pos=$pos&ref=$ref&alt=$alt&genbankid='+document.getElementById('inputString').value;"/>
      <br>
		 
</form>
		
    </div>
END_HTML
    
    print <<END_HTML;
    </body>
    </html>
END_HTML
}

1;
