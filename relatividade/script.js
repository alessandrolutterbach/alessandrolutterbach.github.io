let toggleNavStatus = false;

let getSidebar = document.querySelector('.nav-sidebar');
let getSidebarUl = document.querySelector('.nav-sidebar ul');
let getSidebarLinks = document.querySelector('.nav-sidebar a');

getSidebar.style.width = "0px";
getSidebar.style.backgroundColor = "#ffff";

getSidebarUl.style.visibility = "hidden";

let toggleNav = function() {
    
    if(toggleNavStatus === false) {
        getSidebarUl.style.visibility = "visible";
        getSidebar.style.width = "272px";
        getSidebar.style.backgroundColor = "#1b1b1b";

        let arrayLength = getSidebarLinks.length;
        
        for(let i = 0; i < arrayLength; i++ ) {
            getSidebarLinks[i].style.opacity = "1";
        }

        toggleNavStatus = true;
     } else if(toggleNavStatus === true) {
        getSidebar.style.width = "0px";
        getSidebar.style.backgroundColor = "#ffff";

        let arrayLength = getSidebarLinks.length;
        for(let i = 0; i < arrayLength; i++ ) {
            getSidebarLinks[i].style.opacity = "0";
        }
        
        getSidebarUl.style.visibility = "hidden";
        toggleNavStatus = false;
     }
}