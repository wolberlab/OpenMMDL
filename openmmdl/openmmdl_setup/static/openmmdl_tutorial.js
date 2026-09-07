(function () {
  "use strict";

  var config = window.OpenMMDLTutorial;
  if (!config || !config.steps || !config.steps.length) return;

  // Always begin each page's tutorial at its first step.
  var currentIndex = 0;

  function getTarget(step) {
    if (!step || !step.target) return null;
    try {
      return document.querySelector(step.target);
    } catch (e) {
      return null;
    }
  }

  function clearExpectedInputHint() {
    document.querySelectorAll(".om-tutorial-input-hint, .om-tutorial-option-hint").forEach(function (node) {
      node.remove();
    });
  }

  function showExpectedInputHint(step, target) {
    clearExpectedInputHint();

    if (!step || !target) return;

    var label = step.expectedInputLabel;
    if (!label) return;

    var hint = document.createElement("div");
    hint.className = "om-tutorial-input-hint";
    hint.textContent = label;

    /*
    * Additional-molecule rows: place the hint below the whole companion block
    * (upload + remove + topology code), not inside the flex row.
    */
    if (target.closest) {
      var companionBlock = target.closest(".om-companion-block");
      if (companionBlock && companionBlock.parentNode) {
        companionBlock.parentNode.insertBefore(hint, companionBlock.nextSibling);
        return;
      }
    }

    /*
    * File inputs are wrapped in:
    * .col-md-10 > .om-filedrop
    *
    * Insert the hint directly after .om-filedrop inside that same column.
    */
    if (target.classList && target.classList.contains("om-filedrop")) {
      target.parentNode.insertBefore(hint, target.nextSibling);
      return;
    }

    /*
    * If the actual target is the hidden file input, climb to .om-filedrop first.
    */
    if (target.closest) {
      var fileDrop = target.closest(".om-filedrop");
      if (fileDrop && fileDrop.parentNode) {
        fileDrop.parentNode.insertBefore(hint, fileDrop.nextSibling);
        return;
      }
    }

    /*
    * Radio/checkbox cards (.om-choice-card) wrap the input together with its
    * description, so place the hint below the whole card instead of inside it.
    */
    if (target.closest) {
      var choiceCard = target.closest(".om-choice-card");
      if (choiceCard && choiceCard.parentNode) {
        choiceCard.parentNode.insertBefore(hint, choiceCard.nextSibling);
        return;
      }
    }

    /*
    * Checkbox/radio inputs sit before their visible text inside their own
    * <label>, so place the hint after the whole label (not between the box and
    * its text) when the label isn't part of a .form-group.
    */
    if (target.closest && (target.type === "checkbox" || target.type === "radio")) {
      var wrapLabel = target.closest("label");
      if (wrapLabel && wrapLabel.parentNode && !wrapLabel.closest(".form-group")) {
        wrapLabel.parentNode.insertBefore(hint, wrapLabel.nextSibling);
        return;
      }
    }

    /*
    * Plain buttons (e.g. "+ Add additional molecule"): anchor the hint directly
    * after the button so it sits next to the control, not below the surrounding
    * form-group's help text.
    */
    if (
      target.tagName === "BUTTON" ||
      (target.tagName === "INPUT" && (target.type === "button" || target.type === "submit"))
    ) {
      target.insertAdjacentElement("afterend", hint);
      return;
    }

    /*
    * Normal dropdown/input case.
    */
    if (target.closest) {
      var formGroup = target.closest(".form-group");
      if (formGroup && formGroup.parentNode) {
        formGroup.parentNode.insertBefore(hint, formGroup.nextSibling);
        return;
      }
    }

    target.insertAdjacentElement("afterend", hint);
  }

  function highlightElement(element) {
    clearExpectedOptionHint();
    document.querySelectorAll(".om-tutorial-highlight").forEach(function (node) {
      node.classList.remove("om-tutorial-highlight");
    });
    if (!element) return;
    var visibleElement = element;
    if (element.type === "hidden" || element.offsetParent === null) {
      visibleElement = element.closest(".om-filedrop, .form-group, label, table") || element;
    }
    visibleElement.classList.add("om-tutorial-highlight");
    try {
      visibleElement.scrollIntoView({ behavior: "smooth", block: "center", inline: "nearest" });
    } catch (e) {
      visibleElement.scrollIntoView();
    }
  }

  function clearExpectedOptionHint() {
    document.querySelectorAll(".om-tutorial-option-hint").forEach(function (node) {
      node.remove();
    });
  }

  function showExpectedOptionHint(step, target) {
    clearExpectedOptionHint();

    if (!step || !step.expectedOptionLabel || !target || target.tagName !== "SELECT") return;

    var hint = document.createElement("div");
    hint.className = "om-tutorial-option-hint";
    hint.textContent = "Choose: " + step.expectedOptionLabel;

    target.insertAdjacentElement("afterend", hint);
  }

  function setValue(selector, value) {
    var element = document.querySelector(selector);
    if (!element) return;
    element.value = value;
    triggerChange(element);
  }

  function setChecked(selector, checked) {
    var element = document.querySelector(selector);
    if (!element) return;
    element.checked = checked;
    triggerChange(element);
  }

  function triggerChange(element) {
    if (!element) return;
    element.dispatchEvent(new Event("input", { bubbles: true }));
    element.dispatchEvent(new Event("change", { bubbles: true }));
  }

  // A step may declare a "done" condition; when it holds, the panel shows a
  // check mark. Steps are never advanced automatically: the user clicks Next.
  function stepDone(step) {
    var done = step && step.done;
    if (!done) return false;
    var element;
    try {
      element = document.querySelector(done.selector || step.target);
    } catch (e) {
      return false;
    }
    if (!element) return false;
    if (done.condition === "exists") return true;
    if (done.condition === "checked") return !!element.checked;
    if (done.condition === "hasValue") return !!String(element.value || "").trim();
    if (done.condition === "valueEquals") return String(element.value) === String(done.value);
    return false;
  }

  function updateDoneBadge() {
    var badge = document.querySelector("[data-om-tutorial-done]");
    if (!badge) return;
    badge.hidden = !stepDone(config.steps[currentIndex]);
  }

  function refreshPageUi() {
    if (typeof window.optionChanged === "function") window.optionChanged();
    if (typeof window.waterBoxChanged === "function") window.waterBoxChanged();
  }

  function runAction(action) {
    if (!action) return;

    if (action === "selectPdbPath") {
      setChecked("input[name='type'][value='pdb']", true);
      return;
    }

    if (action === "configurePdbTutorialDefaults") {
      setValue("#smallMoleculeMode", "single");
      setValue("#sdfResname", "UNK");
      refreshPageUi();
      return;
    }

    if (action === "setForceFieldAmber19") {
      setValue("#forcefield", "AMBER19");
      refreshPageUi();
      return;
    }

    if (action === "setWaterModelOpc3") {
      setValue("#waterModel", "OPC3");
      refreshPageUi();
      return;
    }

    if (action === "addTutorialCompanionMolecule") {
      // Add only a single companion row for the tutorial; don't stack duplicates
      // if the button is pressed again.
      var container = document.querySelector("#companionFilesContainer");
      if (container && container.querySelector("[data-companion-row]")) return;
      var addBtn = document.querySelector("#addCompanionMolecule");
      if (addBtn) addBtn.click();
      return;
    }

    if (action === "setLigandForceFieldSmirnoff") {
      setValue("#smallMoleculeForceField", "smirnoff");
      refreshPageUi();
      return;
    }

    if (action === "selectFirstTwoChains") {
      var checkboxes = document.querySelectorAll("#table input[name='include']");
      checkboxes.forEach(function (checkbox, index) {
        checkbox.checked = index < 2;
      });
      setValue("#heterogens", "all");
      refreshPageUi();
      return;
    }

    if (action === "selectAllMissingResidues") {
      var addBoxes = document.querySelectorAll("#table input[name='add']");
      addBoxes.forEach(function (checkbox) {
        checkbox.checked = true;
      });
      return;
    }

    if (action === "setTutorialPh") {
      setChecked("#addHydrogens", true);
      setValue("#phfield", "7.4");
      refreshPageUi();
      return;
    }

    if (action === "enableTutorialWaterBox") {
      setChecked("#addMembrane", false);
      setChecked("#addWater", true);
      setValue("#geomPadding", "1.0");
      setValue("#ionicstrengthfield", "0.15");
      refreshPageUi();
      return;
    }

    if (action === "setSimulationLength") {
      setValue("#sim_length", "10");
      refreshPageUi();
    }
  }

  function buildPanel() {
    var panel = document.createElement("aside");
    panel.className = "om-tutorial-panel";
    panel.setAttribute("aria-live", "polite");
    document.body.appendChild(panel);
    return panel;
  }

  function render(panel) {
    var step = config.steps[currentIndex];
    var target = getTarget(step);
    highlightElement(target);
    showExpectedInputHint(step, target);
    showExpectedOptionHint(step, target);

    var filesLink = config.filesUrl
      ? '<a class="btn btn-default btn-xs" href="' + config.filesUrl + '">Download files</a>'
      : "";
    var docsLink = config.docsUrl
      ? '<a class="btn btn-default btn-xs" href="' + config.docsUrl + '" target="_blank" rel="noopener">Open docs</a>'
      : "";
    var actionButton = step.action
      ? '<button type="button" class="btn btn-primary btn-sm" data-om-tutorial-action="' + step.action + '">' +
        (step.actionLabel || "Apply") +
        "</button>"
      : "";

    // On the final step of the final page, the exit becomes "Finish tutorial".
    var isTutorialEnd =
      config.pageNumber === config.totalPages && currentIndex === config.steps.length - 1;
    var exitLabel = isTutorialEnd ? "Finish tutorial" : "Exit";

    panel.innerHTML =
      '<div class="om-tutorial-kicker">' + config.title + " · page " + config.pageNumber + " of " + config.totalPages + "</div>" +
      '<h3 class="om-tutorial-heading">' + config.pageTitle + "</h3>" +
      '<div class="om-tutorial-progress"><span style="width:' + ((currentIndex + 1) / config.steps.length * 100) + '%"></span></div>' +
      '<div class="om-tutorial-step-count">Step ' + (currentIndex + 1) + " of " + config.steps.length + " on this page</div>" +
      '<h4 class="om-tutorial-step-title">' + step.title +
        '<span class="om-tutorial-done" data-om-tutorial-done hidden>&#10003; Done</span></h4>' +
      '<p class="om-tutorial-body">' + step.body + "</p>" +
      '<div class="om-tutorial-links">' + filesLink + docsLink + "</div>" +
      '<div class="om-tutorial-actions">' +
        '<button type="button" class="btn btn-default btn-sm" data-om-tutorial-prev>Back</button>' +
        '<button type="button" class="btn btn-default btn-sm" data-om-tutorial-next>Next</button>' +
        actionButton +
        '<button type="button" class="btn btn-link btn-sm" data-om-tutorial-restart title="Go back to the first step of this page">Restart page</button>' +
        '<a class="btn btn-link btn-sm" href="' + config.stopUrl + '">' + exitLabel + "</a>" +
      "</div>";

    var prevButton = panel.querySelector("[data-om-tutorial-prev]");
    var nextButton = panel.querySelector("[data-om-tutorial-next]");
    var restartButton = panel.querySelector("[data-om-tutorial-restart]");
    prevButton.disabled = currentIndex === 0;
    nextButton.disabled = currentIndex === config.steps.length - 1;
    restartButton.disabled = currentIndex === 0;

    restartButton.addEventListener("click", function () {
      currentIndex = 0;
      render(panel);
    });

    prevButton.addEventListener("click", function () {
      currentIndex = Math.max(0, currentIndex - 1);
      render(panel);
    });

    nextButton.addEventListener("click", function () {
      currentIndex = Math.min(config.steps.length - 1, currentIndex + 1);
      render(panel);
    });

    var applyButton = panel.querySelector("[data-om-tutorial-action]");
    if (applyButton) {
      applyButton.addEventListener("click", function () {
        runAction(applyButton.getAttribute("data-om-tutorial-action"));
        updateDoneBadge();
      });
    }

    updateDoneBadge();
  }

  function enforceSingleCompanionMolecule() {
    var container = document.querySelector("#companionFilesContainer");
    var addBtn = document.querySelector("#addCompanionMolecule");
    if (!container || !addBtn) return;

    function sync() {
      var row = container.querySelector("[data-companion-row]");
      addBtn.disabled = !!row;
      addBtn.title = row
        ? "The tutorial uses a single additional molecule. Remove it to add a different one."
        : "";
      updateDoneBadge();
    }

    new MutationObserver(sync).observe(container, { childList: true });
    sync();
  }

  document.addEventListener("DOMContentLoaded", function () {
    var panel = buildPanel();
    render(panel);
    enforceSingleCompanionMolecule();
    document.addEventListener("change", updateDoneBadge, true);
    document.addEventListener("input", updateDoneBadge, true);
    document.addEventListener("click", updateDoneBadge, true);
  });
})();
