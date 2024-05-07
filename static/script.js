$(document).ready(function(){
    $("#hide").click(function(){
        $(".cache").hide('slow');
    });
    $("#show").click(function(){
        $(".cache").show('slow');
    });
	
});

jQuery(document).ready(function($){
  // Get current path and find target link
  var path = window.location.pathname.split("/").pop();
  
  // Account for home page with empty path
  if ( path == '' ) {
    path = 'index.php';
  }
      
  var target = $('#bloc-menu li a[href="'+path+'"]');
  // Add active class to target link
  target.addClass('active');
});