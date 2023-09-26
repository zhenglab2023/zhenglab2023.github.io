/***
This class controls intializing and animating on repeat
the cells graphic on the homepage
***/
function CellsAnimation() {
  const numPerRow = [11, 19, 22, 28, 21, 17];
  //const numPerRow = [11, 19, 22, 28, 21, 11];
  const classes = ["cg-cell--red", "cg-cell--blue", "cg-cell--teal", "cg-cell--green", "cg-cell--yellow", "cg-cell--lightblue"];
  var randomCellType = function () { return Math.floor(Math.random()*classes.length) }

  // keeps track of how many cycles we've done
  var turn = 0;
  var itor = 0;

  this.init = function () {
    // Assign a true cell type to each gray cell
    $('.cg-cell').each(function() {
      var type = randomCellType();
      $(this).attr("data-type", type);
      $(this).addClass(classes[type]);
    });

    // This changes the frequency of the animation
    // Commenting this out will make animation not happen and cells stay gray
    setInterval(this.animate, 2500);
    // $('.cg-cell[data-type=2]').addClass('revealed'); //test
  }

  this.animate = function () {
    console.log("Animating");

    if (turn % 2 == 0) {
      // Choose a random cell type to reveal
	var classType = itor;//randomCellType();
      $(`.cg-cell[data-type=${classType}]`).addClass("revealed");
      itor++;
      if(itor >= classes.length) {itor = 0;}
    } else {
      // Return all cells to gray
      $('.cg-cell').removeClass('revealed');
    }
    turn++;
  }
}

$(document).ready(function() {
  var animation = new CellsAnimation();
  animation.init();
});