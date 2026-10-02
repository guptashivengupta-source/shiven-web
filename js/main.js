// menu, punnett square, free throw slider, prime spiral and copy button

const menuBtn = document.querySelector('.menu-toggle');
const nav = document.querySelector('#site-nav');

menuBtn.addEventListener('click', () => {
  const open = menuBtn.getAttribute('aria-expanded') === 'true';
  menuBtn.setAttribute('aria-expanded', !open);
  nav.classList.toggle('is-open', !open);
});

nav.addEventListener('click', (e) => {
  if (e.target.tagName === 'A') {
    menuBtn.setAttribute('aria-expanded', false);
    nav.classList.remove('is-open');
  }
});


const parent1 = document.querySelector('#parent-1');
const parent2 = document.querySelector('#parent-2');
const grid = document.querySelector('#punnett-grid');
const result = document.querySelector('#punnett-result');

function updatePunnett() {
  const a = parent1.value.split('');
  const b = parent2.value.split('');
  const count = { PP: 0, Pp: 0, pp: 0 };
  let rows = '';

  for (const x of a) {
    let cells = '';
    for (const y of b) {
      const kid = [x, y].sort().join('');
      count[kid]++;
      cells += `<td class="${kid.includes('P') ? 'is-purple' : 'is-white'}">${kid}</td>`;
    }
    rows += `<tr><th scope="row">${x}</th>${cells}</tr>`;
  }

  grid.innerHTML = `
    <caption class="visually-hidden">Offspring of ${parent1.value} x ${parent2.value}</caption>
    <thead><tr><td></td>${b.map(y => `<th scope="col">${y}</th>`).join('')}</tr></thead>
    <tbody>${rows}</tbody>`;

  const types = Object.keys(count)
    .filter(k => count[k] > 0)
    .map(k => `${count[k]} <span class="alleles">${k}</span>`);
  const purple = (count.PP + count.Pp) * 25;

  result.innerHTML = `
    <p><strong>Genotypes</strong>${types.join(' : ')}</p>
    <p><strong>Flowers</strong>${purple}% purple, ${100 - purple}% white</p>`;
}

parent1.addEventListener('change', updatePunnett);
parent2.addEventListener('change', updatePunnett);
document.querySelector('#punnett-form').addEventListener('submit', e => e.preventDefault());


const slider = document.querySelector('#ft-rate');

slider.addEventListener('input', () => {
  const p = slider.value / 100;
  const points = p + p * p;
  document.querySelector('#ft-rate-value').textContent = slider.value + '%';
  document.querySelector('#ft-worth').textContent = points.toFixed(2) + ' points';
  document.querySelector('#ft-verdict').textContent = points > 1 ? 'more than a sure point' : 'less than a sure point';
});


const canvas = document.querySelector('#prime-spiral');

function drawSpiral() {
  const size = 201;
  const total = size * size;

  const prime = new Array(total + 1).fill(true);
  prime[0] = prime[1] = false;
  for (let i = 2; i * i <= total; i++) {
    if (!prime[i]) continue;
    for (let j = i * i; j <= total; j += i) prime[j] = false;
  }

  const dpr = window.devicePixelRatio || 1;
  const space = Math.min(canvas.parentElement.clientWidth || 220, 220);
  const cell = Math.max(1, Math.floor(space * dpr / size));
  canvas.width = canvas.height = cell * size;
  canvas.style.width = canvas.style.height = (cell * size / dpr) + 'px';

  const ctx = canvas.getContext('2d');
  ctx.fillStyle = getComputedStyle(canvas).color;

  const dirs = [[1, 0], [0, -1], [-1, 0], [0, 1]];
  let x = (size - 1) / 2;
  let y = x;
  let n = 1;
  let step = 1;
  let d = 0;

  while (n < total) {
    for (let k = 0; k < 2; k++) {
      for (let i = 0; i < step && n < total; i++) {
        x += dirs[d][0];
        y += dirs[d][1];
        n++;
        if (prime[n]) ctx.fillRect(x * cell, y * cell, cell, cell);
      }
      d = (d + 1) % 4;
    }
    step++;
  }
}

drawSpiral();
window.addEventListener('resize', drawSpiral);
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', drawSpiral);
new MutationObserver(drawSpiral).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });


const copyBtn = document.querySelector('.copy-btn');

copyBtn.addEventListener('click', async () => {
  try {
    await navigator.clipboard.writeText(copyBtn.dataset.copy);
    copyBtn.textContent = 'Copied';
  } catch (err) {
    getSelection().selectAllChildren(document.querySelector('#email'));
    copyBtn.textContent = 'Press Ctrl+C';
  }
  setTimeout(() => copyBtn.textContent = 'Copy', 2000);
});
