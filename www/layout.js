// Fit both plots and the settings body to the actual viewport.
(() => {
  function initialize() {
    const cards = [...document.querySelectorAll('.workbench > .map, .workbench > .detail-map')];
    const main = document.querySelector('main');
    const controls = document.querySelector('.controls');
    const body = document.querySelector('.controls-body');
    if (!cards.length || !main) return;
    let frame = 0;
    function fit() {
      frame = 0;
      const bottomSpace = Math.max(18, parseFloat(getComputedStyle(main).paddingBottom) || 0);
      const headingHeight = Math.ceil(Math.max(...cards.map(card => card.querySelector('.map-heading-content').getBoundingClientRect().height)));
      cards.forEach(card => {
        const heading = card.querySelector('.map-heading');
        if (heading.style.minHeight !== `${headingHeight}px`) heading.style.minHeight = `${headingHeight}px`;
      });
      const top = cards[0].getBoundingClientRect().top + scrollY;
      const chrome = Math.max(...cards.map(card => card.getBoundingClientRect().height - card.querySelector('.html-widget-output').getBoundingClientRect().height));
      const height = `${Math.max(300, Math.floor(innerHeight - top - chrome - bottomSpace))}px`;
      cards.forEach(card => { if (card.style.getPropertyValue('--umap-fit-height') !== height) card.style.setProperty('--umap-fit-height', height); });
      if (controls && body) {
        const r = controls.getBoundingClientRect();
        const limit = `${Math.max(200, Math.floor(innerHeight - r.top - scrollY - (r.height - body.getBoundingClientRect().height) - bottomSpace))}px`;
        if (controls.style.getPropertyValue('--controls-fit-height') !== limit) controls.style.setProperty('--controls-fit-height', limit);
      }
    }
    function schedule() { if (!frame) frame = requestAnimationFrame(fit); }
    const observer = new ResizeObserver(schedule);
    [...cards, controls, document.querySelector('.top'), document.querySelector('.dataset-status')].filter(Boolean).forEach(e => observer.observe(e));
    addEventListener('resize', schedule);
    document.addEventListener('visibilitychange', schedule);
    if (window.jQuery) window.jQuery(document).on('shiny:connected shiny:value shown.bs.tab', schedule);
    if (document.fonts) document.fonts.ready.then(schedule);
    document.addEventListener('click', event => {
      const button = event.target.closest('button[data-annotation]');
      if (button && window.Shiny) Shiny.setInputValue('annotation_pick', button.dataset.annotation, {priority: 'event'});
    });
    schedule();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initialize); else initialize();
})();
