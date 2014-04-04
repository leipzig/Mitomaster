$(function(){
      $('.tooltip-species').tooltip()
      
      $('.popover-test').popover()

      // popover demo
      $("a[rel=popover]")
        .popover()
        .click(function(e) {
          e.preventDefault()
        })
    
  
        //all species in all cats
        $("A[href='#select_everything']").click( function() {
    		//without the change() event the related ajax doesn't occur
    		//however it sends the whole form for every checkbox
              $("#all_species INPUT[type='checkbox']").attr('checked', true);
              //.change();
              send_species_form();
              return false;
          });
          
          //clear all except rCRS
          $("A[href='#select_crs']").click( function() {
      		//without the change() event the related ajax doesn't occur
      		//however it sends the whole form for every checkbox
                $("#all_species INPUT[type='checkbox']").attr('checked', false);
                $("#Homo_sapiens_\\(rCRS\\)").attr('checked', true);
                //.change();
                send_species_form();
                return false;
            });
       
       //primates only
       $("A[href='#select_primates']").click( function() {
     		//without the change() event the related ajax doesn't occur
     		//however it sends the whole form for every checkbox
               $("#Primates INPUT[type='checkbox']").attr('checked', true);
               $("#Homo_sapiens_\\(rCRS\\)").attr('checked', true);
               //.change();
               send_species_form();
               return false;
           });       
           
      // Select all in a cateogry
      $("A[href='#select_all']").click( function() {
		//without the change() event the related ajax doesn't occur
		//however it sends the whole form for every checkbox
          $("#" + $(this).attr('rel') + " INPUT[type='checkbox']").attr('checked', true)
          //.change();
          send_species_form();
          return false;
      });

      // Select none
      $("A[href='#select_none']").click( function() {
          $("#" + $(this).attr('rel') + " INPUT[type='checkbox']").attr('checked', false)
          //.change();
          send_species_form();
          return false;
      });
      
      $('.genbank_link').live('click', function() {
          $("#gb_accession").val($(this).html());
      });

      var snv_list = new Array();
      snv_list['snv_tabular']="sample	pos	ref	var\nSample1	73	A	G\nSample1	146	T	C\nSample2	263	A	:\nSample2	709	:	G"
      snv_list['snv_compact']="Sample1	A73G\nSample1	T146C\nSample2	A263d\nSample2	709.1G"
      
      $('.snv_link').live('click', function() {
          // http://stackoverflow.com/questions/415602/set-value-of-textarea-in-jquery
          $("textarea#snv_list").val(snv_list[$(this).attr("id")]);
      });
      
      send_species_form = function() {
          $("#success").spin()
           var formData = $("#species_checkbox").serialize();
  		//console.log(formData);
           $.ajax({
              url: "cons_checkboxes.cgi",
              type: "POST",
              data: formData,
              cache: false,
              success: function(data){
                  console.log(data);
                  $("#success").spin(false);
                },
              dataType: "json"
        });
    }
      
      $(".species_check").change(send_species_form);
      
      //teh spin.js jquery plugin
      $.fn.spin = function(opts) {
        this.each(function() {
          var $this = $(this),
              data = $this.data();

          if (data.spinner) {
            data.spinner.stop();
            delete data.spinner;
          }
          if (opts !== false) {
            data.spinner = new Spinner($.extend({color: $this.css('color')}, opts)).spin(this);
          }
        });
        return this;
      };
});

