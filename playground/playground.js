(function () {
  "use strict";

  var input = document.getElementById("formula");
  var output = document.getElementById("output");
  var prompt = document.querySelector(".prompt");
  var modes = document.querySelectorAll("input[name=mode]");
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

  function mode() {
    var checked = document.querySelector("input[name=mode]:checked");
    return checked ? checked.value : "formula";
  }

  // One input is read at a time; while it runs, only the latest edit waits.
  function request(message) {
    if (worker === null) {
      return;
    }
    if (busy) {
      queued = message;
      return;
    }
    busy = true;
    output.setAttribute("aria-busy", "true");
    worker.postMessage(message);
  }

  function update() {
    var text = input.value.trim();
    prompt.textContent = mode() === "smiles" ? "SMILES" : "Formula";
    if (text === "") {
      queued = null;
      output.removeAttribute("aria-busy");
      say(mode() === "smiles" ? "Type a SMILES string to draw it." : "Type a formula to read it.");
      return;
    }
    request({ input: text, mode: mode() });
  }

  function choose(value, text) {
    Array.prototype.forEach.call(modes, function (radio) {
      radio.checked = radio.value === value;
    });
    input.value = text;
    update();
    input.focus();
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
    if (reply.input !== input.value.trim() || reply.mode !== mode()) {
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
  Array.prototype.forEach.call(modes, function (radio) {
    radio.addEventListener("change", update);
  });
  Array.prototype.forEach.call(document.querySelectorAll(".examples button"), function (button) {
    button.addEventListener("click", function () {
      if (button.dataset.smiles) {
        choose("smiles", button.dataset.smiles);
      } else {
        choose("formula", button.textContent);
      }
    });
  });
  if (worker !== null) {
    update();
  }
})();
