function init() {
    if (localStorage.getItem("mode") == "light") {
        light();
    }
    var element = document.body;
    element.classList.add("selected1");
}

function light() {
    var element = document.body;
    element.classList.add("dark-mode");
    document.getElementById("potencialluz").src="img/potencialluz2.png";
}

function change() {
    var element = document.body;
    element.classList.toggle("dark-mode");
    if (localStorage.getItem("mode") == "dark") {
        localStorage.setItem("mode", "light");
        document.getElementById("potencialluz").src="img/potencialluz2.png";
    } else {
        localStorage.setItem("mode", "dark");
        document.getElementById("potencialluz").src="img/potencialluz.png";
    }
}

function show(shown, hidden, hide, select) {
    document.getElementById(shown).style.display='block';
    document.getElementById(hidden).style.display='none';
    document.getElementById(hide).style.display='none';
    var element = document.body;
    if (select == "um") {
        element.classList.add("selected1");
        element.classList.remove("selected2");
        element.classList.remove("selected3");
    } else if (select == "dois") {
        element.classList.add("selected2");
        element.classList.remove("selected1");
        element.classList.remove("selected3");
    } else if (select == "tres") {
        element.classList.add("selected3");
        element.classList.remove("selected2");
        element.classList.remove("selected1");
    }
    return false
}

function getValues() {
    m = document.getElementById("rangevalue").value;
    e = document.getElementById("number").value;
    o = document.getElementById("rangevalue2").value;
    mode = localStorage.getItem("mode");
    D1 = document.getElementById("rangevalue3").value;
    D2 = document.getElementById("number2").value;
}

setInterval(getValues,100);

function limpaDiv(show){
    document.getElementById(show).style.display='block';
    var graph = document.getElementById("graph");
    graph.firstChild.remove();
    graph.firstChild.remove();
    return false
}

function limpaDiv2() {
    var graph = document.getElementById("graph2");
    graph.firstChild.remove();
}

function limpaDiv3() {
    var graph = document.getElementById("graph3");
    graph.firstChild.remove();
}