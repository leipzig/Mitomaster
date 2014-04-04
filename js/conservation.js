$.getJSON( "http://wncg.ece.utexas.edu/intranet/admin/list.php?was=<?php echo $lookup;?>", null, function (json) {
    $("#list").dataTable( {
        "bProcessing": true,
        "bPaginate": true,
        "sPaginationType": "full_numbers",
        "aaData": json.aaData,
        "aoColumns": json.aoColumns
    } );
} );