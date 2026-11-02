(function () {
  "use strict";

  var chars = Array.prototype.slice.call(document.querySelectorAll(".formula [data-offset]"));
  var parts = Array.prototype.slice.call(document.querySelectorAll("svg [data-span]"));
  var selected = null;

  function spanOf(el) {
    if (el.hasAttribute("data-span")) {
      return el.getAttribute("data-span").split("-").map(Number);
    }
    return [Number(el.getAttribute("data-offset")), Number(el.getAttribute("data-end"))];
  }

  function overlaps(a, b) {
    return a[0] < b[1] && b[0] < a[1];
  }

  // A part lights up its own characters only; characters light up every part they produced.
  function show(span, part) {
    chars.forEach(function (c) {
      c.classList.toggle("hl", span !== null && overlaps(spanOf(c), span));
    });
    parts.forEach(function (p) {
      var lit = span !== null && (part ? p === part : overlaps(spanOf(p), span));
      p.classList.toggle("hl", lit);
    });
  }

  function traced(event) {
    var el = event.target;
    return el && el.closest ? el.closest("svg [data-span], .formula [data-offset]") : null;
  }

  function enter(event) {
    var el = traced(event);
    if (el) {
      show(spanOf(el), el.hasAttribute("data-span") ? el : null);
    }
  }

  function leave(event) {
    if (traced(event)) {
      show(selected, null);
    }
  }

  function selection() {
    var current = document.getSelection();
    if (!current || current.rangeCount === 0 || current.isCollapsed) {
      return null;
    }
    var range = current.getRangeAt(0);
    var picked = chars.filter(function (c) {
      return range.intersectsNode(c);
    });
    if (picked.length === 0) {
      return null;
    }
    return [spanOf(picked[0])[0], spanOf(picked[picked.length - 1])[1]];
  }

  chars.concat(parts).forEach(function (el) {
    el.setAttribute("tabindex", "0");
  });
  document.querySelectorAll("svg.moiety-structure").forEach(function (svg) {
    svg.setAttribute("role", "group");
  });
  document.addEventListener("mouseover", enter);
  document.addEventListener("mouseout", leave);
  document.addEventListener("focusin", enter);
  document.addEventListener("focusout", leave);
  document.addEventListener("selectionchange", function () {
    var next = selection();
    if (next !== null || selected !== null) {
      selected = next;
      show(selected, null);
    }
  });
})();
