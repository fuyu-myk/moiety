(function () {
  "use strict";

  var input = document.getElementById("formula");
  var output = document.getElementById("output");
  var worker = null;
  var busy = false;
  var queued = null;

  function say(text) {
    var summary = document.createElement("p");
    summary.className = "summary";
    summary.textContent = text;
    output.replaceChildren(summary);
  }

  function stop(reason) {
    if (worker !== null) {
      worker.terminate();
      worker = null;
    }
    busy = false;
    queued = null;
    output.removeAttribute("aria-busy");
    say("The playground could not start: " + reason);
  }

  // One formula is read at a time; while it runs, only the latest edit waits.
  function request(formula) {
    if (worker === null) {
      return;
    }
    if (busy) {
      queued = formula;
      return;
    }
    busy = true;
    output.setAttribute("aria-busy", "true");
    worker.postMessage(formula);
  }

  function update() {
    var formula = input.value.trim();
    if (formula === "") {
      queued = null;
      output.removeAttribute("aria-busy");
      say("Type a formula to read it.");
      return;
    }
    request(formula);
  }

  function receive(event) {
    var reply = event.data;
    var next = queued;
    busy = false;
    queued = null;
    if (next !== null) {
      request(next);
    } else {
      output.removeAttribute("aria-busy");
    }
    if (reply.formula !== input.value.trim()) {
      return;
    }
    if (reply.error) {
      say("moiety could not run: " + reply.error);
    } else {
      output.innerHTML = reply.html;
    }
  }

  try {
    worker = new Worker("worker.js");
    worker.onmessage = receive;
    worker.onerror = function (event) {
      stop(event.message || "worker.js did not load");
    };
  } catch (error) {
    stop(String(error));
  }

  input.addEventListener("input", update);
  Array.prototype.forEach.call(document.querySelectorAll(".examples button"), function (button) {
    button.addEventListener("click", function () {
      input.value = button.textContent;
      update();
      input.focus();
    });
  });
  if (worker !== null) {
    update();
  }
})();
