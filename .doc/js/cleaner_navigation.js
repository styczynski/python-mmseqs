(function(){
    $(document).ready(function(){
        $(".sidebar-tree a").each(function() {
           const text = $(this).text();
           if (text.indexOf("biosnake.") > -1) {
              const tokens = text.split('.');
              $(this).text(tokens[tokens.length-1]);
           }
        });
    });
})();