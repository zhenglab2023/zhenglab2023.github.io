$(document).ready(function() {
  const link = "mail" + "to:";
  const domain = "@" + "stanford.edu";

  $(".email").each(function() {
    var name = $(this).data("name");
    var address = name + domain;
    $(this).attr("href", link + address);
    $(this).html(address);
  });

  $(".collapser").click(function() {
    $(this).siblings(".collapse").collapse("toggle");
    return false;
  });
});