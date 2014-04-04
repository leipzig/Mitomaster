$(function(){
    
    $('#gb_accession').typeahead({
        source: mitoids
    });

    //this prevents the autocomplete from timing out for slow-clickers
    //https://github.com/twitter/bootstrap/issues/2715

    
});