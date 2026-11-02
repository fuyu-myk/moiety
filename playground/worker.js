"use strict";

// Reads formulas off the main thread, so a slow formula never blocks typing.
// Each message is a formula; each reply is { formula, html } or { formula, error }.

// Wasm memory never shrinks, so an instance that grew past this is replaced.
var MEMORY_LIMIT = 64 * 1024 * 1024;

var encoder = new TextEncoder();
var decoder = new TextDecoder();
var compiled = null;
var instance = null;

// A failed load is forgotten, so the next formula tries again.
function load() {
  if (compiled === null) {
    compiled = fetch("moiety.wasm")
      .then(function (response) {
        if (!response.ok) {
          throw new Error("moiety.wasm: HTTP " + response.status);
        }
        return response.arrayBuffer();
      })
      .then(function (bytes) {
        return WebAssembly.compile(bytes);
      });
    compiled.catch(function () {
      compiled = null;
    });
  }
  return compiled;
}

// Each call comes before the view it needs: calls can grow memory, which detaches older views.
function render(exports, formula) {
  var bytes = encoder.encode(formula);
  var at = exports.moiety_input(bytes.length);
  new Uint8Array(exports.memory.buffer, at, bytes.length).set(bytes);
  var length = exports.moiety_render();
  at = exports.moiety_output();
  return decoder.decode(new Uint8Array(exports.memory.buffer, at, length));
}

self.onmessage = function (event) {
  var formula = event.data;
  load()
    .then(function (module) {
      return instance || WebAssembly.instantiate(module, {});
    })
    .then(function (ready) {
      instance = ready;
      var html = render(instance.exports, formula);
      if (instance.exports.memory.buffer.byteLength > MEMORY_LIMIT) {
        instance = null;
      }
      self.postMessage({ formula: formula, html: html });
    })
    .catch(function (error) {
      // A trap can leave the module half way through an update, so the next formula gets a fresh one.
      instance = null;
      self.postMessage({ formula: formula, error: String(error) });
    });
};
