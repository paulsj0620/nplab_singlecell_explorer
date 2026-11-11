// Keep Shiny/Bootstrap's existing tab activation path for clicks and drag release.
(() => {
  function initialize() {
    const nav = document.querySelector('#stage');
    if (!nav || nav.classList.contains('nav-enhanced')) return;
    const indicator = document.createElement('li');
    indicator.className = 'nav-drag-indicator';
    indicator.hidden = true;
    const selection = document.createElement('li');
    selection.className = 'nav-selection-indicator';
    for (const node of [indicator, selection]) {
      node.setAttribute('aria-hidden', 'true');
      node.setAttribute('role', 'presentation');
      nav.append(node);
    }
    let ready = false, gesture = null, suppressClick = false, timer;
    const tabs = () => [...nav.querySelectorAll('a[data-value]')].filter(a =>
      a.getClientRects().length && !a.closest('[hidden], .disabled') && a.getAttribute('aria-disabled') !== 'true');
    function syncSelection() {
      const target = tabs().find(a => a.parentElement.classList.contains('active'));
      if (!target) { selection.hidden = true; ready = false; return; }
      const bounds = nav.getBoundingClientRect(), rect = target.getBoundingClientRect();
      selection.hidden = false;
      selection.style.transition = ready ? '' : 'none';
      Object.assign(selection.style, {
        left: `${rect.left - bounds.left - nav.clientLeft}px`,
        top: `${rect.top - bounds.top - nav.clientTop}px`, width: `${rect.width}px`, height: `${rect.height}px`
      });
      if (!ready) { selection.getBoundingClientRect(); ready = true; selection.style.transition = ''; }
    }
    nav.classList.add('nav-enhanced');
    nav.addEventListener('dragstart', e => e.preventDefault());
    new MutationObserver(records => {
      if (records.some(r => r.target !== selection && r.target !== indicator)) syncSelection();
    }).observe(nav, {subtree: true, attributes: true, attributeFilter: ['class','hidden','style','aria-disabled'], childList: true, characterData: true});
    new ResizeObserver(syncSelection).observe(nav);
    document.addEventListener('languagechange', syncSelection);
    if (window.jQuery) window.jQuery(nav).on('shown.bs.tab', syncSelection);
    syncSelection();
    function preview(x, y) {
      const bounds = nav.getBoundingClientRect(), buttons = tabs();
      const inside = x >= bounds.left && x <= bounds.right && y >= bounds.top - 16 && y <= bounds.bottom + 16;
      // Distance to rectangles also handles navigation wrapped onto multiple rows.
      const distance = a => { const r = a.getBoundingClientRect(); return Math.hypot(Math.max(r.left-x,0,x-r.right), Math.max(r.top-y,0,y-r.bottom)); };
      const target = inside && buttons.length ? buttons.reduce((best,a) => distance(a) < distance(best) ? a : best) : null;
      gesture.target = target;
      indicator.hidden = !target;
      if (!target) return;
      const rect = target.getBoundingClientRect();
      Object.assign(indicator.style, {width: `${rect.width}px`, height: `${rect.height}px`,
        top: `${rect.top-bounds.top-nav.clientTop}px`,
        left: `${Math.max(4,Math.min(nav.clientWidth-rect.width-4,x-bounds.left-nav.clientLeft-rect.width/2))}px`});
    }
    function finish(commit) {
      if (!gesture) return;
      const previous = gesture; gesture = null;
      nav.classList.remove('is-dragging'); indicator.hidden = true;
      if (nav.hasPointerCapture(previous.id)) nav.releasePointerCapture(previous.id);
      if (previous.dragging) {
        suppressClick = true; clearTimeout(timer); timer = setTimeout(() => suppressClick = false, 0);
        if (commit && previous.target) { previous.target.focus({preventScroll:true}); previous.target.click(); }
      }
      syncSelection();
    }
    nav.addEventListener('pointerdown', e => {
      if (!e.isPrimary || e.button !== 0 || !tabs().includes(e.target.closest('a[data-value]'))) return;
      gesture = {id:e.pointerId, x:e.clientX, y:e.clientY, dragging:false, target:null};
    });
    nav.addEventListener('pointermove', e => {
      if (!gesture || e.pointerId !== gesture.id) return;
      if (!gesture.dragging && Math.hypot(e.clientX-gesture.x,e.clientY-gesture.y) < 6) return;
      gesture.dragging = true; nav.setPointerCapture(e.pointerId); nav.classList.add('is-dragging');
      preview(e.clientX,e.clientY); e.preventDefault();
    });
    nav.addEventListener('pointerup', e => { if (gesture?.id !== e.pointerId) return; if (gesture.dragging) preview(e.clientX,e.clientY); finish(true); });
    nav.addEventListener('pointercancel', () => finish(false));
    nav.addEventListener('lostpointercapture', () => finish(false));
    nav.addEventListener('click', e => { if (suppressClick && e.isTrusted && e.detail > 0) { e.preventDefault(); e.stopImmediatePropagation(); } }, true);
    document.addEventListener('keydown', e => { if (e.key === 'Escape' && gesture) { finish(false); e.preventDefault(); } });
    window.addEventListener('blur', () => finish(false));
    window.addEventListener('resize', () => { finish(false); syncSelection(); });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initialize); else initialize();
})();
