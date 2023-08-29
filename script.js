function colorPicker() {
    let tema = ['azul', 'laranja', 'verde', 'roxo', 'amarelo', 'rosa', 'vermelho'];
    let fundos = ['#C9FFFF', '#FFDCC9', '#C9FFD5', '#C9D1FF', '#FFF3C9', '#FFC9FA', '#FFC9C9'];
    let cor = ['#1A5E5E', '#5E271A', '#1A5E25', '#4C1A5E', '#5E571A', '#5E1A57', '#5E1A1A'];

    let min = 0;
    let max = 6;

    var i = Math.floor(Math.random() * (max - min + 1) ) + min;

    var r = document.querySelector(':root');
    r.style.setProperty('--back', fundos[i]);
    r.style.setProperty('--color', cor[i]);

    var b = document.querySelector('body');
    b.classList.add(tema[i]);
}