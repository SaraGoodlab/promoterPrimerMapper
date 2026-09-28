(function () {
  "use strict";
  function updateActions(busy) {
    document.querySelectorAll(".ppm-action").forEach(function (button) {
      button.disabled = busy || button.dataset.ppmEnabled !== "true";
    });
  }
  $(document).on("shiny:connected", function () {
    Shiny.addCustomMessageHandler("ppm-operation", function (message) {
      document.documentElement.dataset.ppmBusy = String(message.active);
      updateActions(message.active);
      var status = document.getElementById("operation_status");
      if (status) status.textContent = message.active ? message.label + "…" : "";
    });
    Shiny.addCustomMessageHandler("ppm-clear-brush", function (message) {
      Shiny.resetBrush("cpg_brush");
    });
    Shiny.addCustomMessageHandler("ppm-validation", function (message) {
      var panel = document.getElementById("advanced_parameters");
      if (panel) panel.open = true;
    });
    Shiny.addCustomMessageHandler("ppm-controls", function (state) {
      Object.keys(state).forEach(function (id) {
        var button = document.getElementById(id);
        if (button) button.dataset.ppmEnabled = String(state[id]);
      });
      updateActions(document.documentElement.dataset.ppmBusy === "true");
    });
  });
  document.addEventListener("click", function (event) {
    var pair = event.target.closest("[data-pair-rank]");
    if (pair) Shiny.setInputValue("pair_row", Number(pair.dataset.pairRank), {priority: "event"});
    var button = event.target.closest("[data-copy-oligo]");
    if (!button) return;
    var status = button.parentElement.querySelector("[role=status]");
    if (!navigator.clipboard) {
      status.textContent = "Select the sequence to copy it.";
      return;
    }
    navigator.clipboard.writeText(button.dataset.copyOligo).then(function () {
      status.textContent = "Copied.";
    }, function () {
      status.textContent = "Copy unavailable. Select the sequence to copy it.";
    });
  });
}());
