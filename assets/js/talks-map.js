// Map of talks and conferences on the Notes page (_includes/talks-map.html).
// Reads _data/talks.yml, inlined as JSON, and draws one pin per city.
(function() {
  var container = document.getElementById('talks-map');
  var dataNode = document.getElementById('talks-map-data');
  if (!container || !dataNode) return;
  var card = container.closest('.talks-map-card');

  function hideCard() {
    if (card) card.hidden = true;
  }

  // Leaflet comes from a CDN; if it was blocked or offline, drop the empty card.
  if (!window.L) {
    hideCard();
    return;
  }

  var data;
  try {
    data = JSON.parse(dataNode.textContent);
  } catch (error) {
    hideCard();
    return;
  }

  var MONTHS = { Jan: 0, Feb: 1, Mar: 2, Apr: 3, May: 4, Jun: 5, Jul: 6, Aug: 7, Sep: 8, Oct: 9, Nov: 10, Dec: 11 };

  // "Oct. 2026" -> a number that sorts by month.
  function dateKey(date) {
    var match = /([A-Z][a-z]{2})\.?\s+(\d{4})/.exec(date || '');
    return match ? Number(match[2]) * 12 + (MONTHS[match[1]] || 0) : 0;
  }

  function escapeHtml(value) {
    return String(value == null ? '' : value).replace(/[&<>"']/g, function(ch) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[ch];
    });
  }

  // Flatten talks and attended conferences into one list of events.
  var events = [];
  (data.papers || []).forEach(function(paper) {
    (paper.talks || []).forEach(function(talk) {
      events.push(Object.assign({ paper: paper.title }, talk));
    });
  });
  (data.attended || []).forEach(function(conference) {
    events.push(Object.assign({ attended: true }, conference));
  });

  // Events in the same city share a pin.
  var places = [];
  var placeByCity = {};
  events.forEach(function(event) {
    if (typeof event.lat !== 'number' || typeof event.lng !== 'number') return;
    var place = placeByCity[event.city];
    if (!place) {
      place = placeByCity[event.city] = { city: event.city, lat: event.lat, lng: event.lng, events: [] };
      places.push(place);
    }
    place.events.push(event);
  });
  places.forEach(function(place) {
    place.events.sort(function(a, b) { return dateKey(b.date) - dateKey(a.date); });
  });

  function pinHtml(place) {
    var count = place.events.length;
    var className = 'talks-map-pin';
    // Hollow until something there has happened; grey when no talk was given there.
    if (place.events.every(function(event) { return event.scheduled; })) className += ' talks-map-pin--scheduled';
    if (place.events.every(function(event) { return event.attended; })) className += ' talks-map-pin--attended';
    if (count > 1) className += ' talks-map-pin--count';
    return '<span class="' + className + '">' + (count > 1 ? count : '') + '</span>';
  }

  var CHEVRON_LEFT = '<svg viewBox="0 0 16 16" width="14" height="14" aria-hidden="true" focusable="false">' +
    '<path d="M10 3.5 5.5 8l4.5 4.5" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></svg>';
  var CHEVRON_RIGHT = '<svg viewBox="0 0 16 16" width="14" height="14" aria-hidden="true" focusable="false">' +
    '<path d="M6 3.5 10.5 8 6 12.5" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></svg>';

  function popupHtml(place) {
    var html = '<p class="talks-map-popup__city">' + escapeHtml(place.city) + '</p>';
    place.events.forEach(function(event) {
      html += '<div class="talks-map-popup__event">';
      if (event.photo) {
        html += '<img class="talks-map-popup__photo" src="' + escapeHtml(event.photo) + '" alt="' +
          escapeHtml(event.photo_alt || event.venue) + '" loading="lazy" decoding="async">';
      }
      html += '<p class="talks-map-popup__venue">' + escapeHtml(event.venue) + '</p>';
      html += '<p class="talks-map-popup__date">' + escapeHtml(event.date) + (event.scheduled ? ' (scheduled)' : '') +
        (event.attended ? ' · Attended' : '') + '</p>';
      if (event.paper) {
        html += '<p class="talks-map-popup__paper"><em>' + escapeHtml(event.paper) + '</em></p>';
      }
      html += '</div>';
    });
    if (place.events.length > 1) {
      html += '<div class="talks-map-popup__pager">' +
        '<button type="button" class="talks-map-popup__step" data-step="-1" aria-label="Previous event">' + CHEVRON_LEFT + '</button>' +
        '<span class="talks-map-popup__count" aria-live="polite"></span>' +
        '<button type="button" class="talks-map-popup__step" data-step="1" aria-label="Next event">' + CHEVRON_RIGHT + '</button>' +
        '</div>';
    }
    return html;
  }

  // A popup shows one event at a time; a city with several events gets a pager
  // (buttons or the left/right arrow keys). It opens on the newest event each time.
  function buildPopup(place, getPopup) {
    var root = document.createElement('div');
    root.className = 'talks-map-popup';
    root.innerHTML = popupHtml(place);
    var slides = root.querySelectorAll('.talks-map-popup__event');
    if (slides.length < 2) return { element: root, reset: function() {} };

    var prev = root.querySelector('[data-step="-1"]');
    var next = root.querySelector('[data-step="1"]');
    var count = root.querySelector('.talks-map-popup__count');
    var current = 0;

    function show(index) {
      current = Math.max(0, Math.min(slides.length - 1, index));
      Array.prototype.forEach.call(slides, function(slide, i) { slide.hidden = i !== current; });
      prev.disabled = current === 0;
      next.disabled = current === slides.length - 1;
      count.textContent = (current + 1) + ' / ' + slides.length;
    }

    function step(delta) {
      var focused = document.activeElement;
      show(current + delta);
      // Re-measure the popup: Leaflet re-attaches the content, which drops focus.
      var popup = getPopup();
      if (popup && popup.isOpen()) popup.update();
      if (focused && focused.disabled) focused = delta < 0 ? next : prev;
      if (focused && root.contains(focused)) focused.focus({ preventScroll: true });
    }

    prev.addEventListener('click', function() { step(-1); });
    next.addEventListener('click', function() { step(1); });
    root.addEventListener('keydown', function(event) {
      if (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight') return;
      // Keep the arrow keys from also panning the map.
      event.preventDefault();
      event.stopPropagation();
      step(event.key === 'ArrowLeft' ? -1 : 1);
    });

    show(0);
    return { element: root, reset: function() { show(0); } };
  }

  var reduceMotion = !!(window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches);
  // Phones, and tablets or touchscreens whose browser reports itself as desktop (iPad Safari).
  var touchFirst = L.Browser.mobile || !!(window.matchMedia && window.matchMedia('(pointer: coarse)').matches);

  var map = L.map(container, {
    // Quarter steps let a phone-width map zoom out far enough to show every pin.
    minZoom: 0,
    // Pins mark cities, not venues, so stop at city scale.
    maxZoom: 11,
    zoomSnap: 0.25,
    worldCopyJump: true,
    // The map only takes the scroll wheel (and, on phones, one-finger drags) once
    // someone has clicked or tapped into it, so scrolling past it never gets stuck.
    scrollWheelZoom: false,
    dragging: !touchFirst,
    // Phones pinch to zoom; on a narrow map the buttons would sit on top of popups.
    zoomControl: !L.Browser.mobile,
    zoomAnimation: !reduceMotion,
    fadeAnimation: !reduceMotion,
    markerZoomAnimation: !reduceMotion
  });
  map.attributionControl.setPrefix('<a href="https://leafletjs.com">Leaflet</a>');

  // CARTO basemaps (English labels): Positron in light mode, Dark Matter in dark mode.
  // The key is carto_basemaps_key in _config.yml; {r} asks for @2x tiles on retina screens.
  var root = document.documentElement;

  function tileStyle() {
    return root.getAttribute('data-theme') === 'dark' ? 'dark_all' : 'light_all';
  }

  var tileLayers = {};

  function tileLayer(style) {
    if (!tileLayers[style]) {
      tileLayers[style] = L.tileLayer('https://basemaps.cartocdn.com/rastertiles/{style}/{z}/{x}/{y}{r}.png?key={key}', {
        style: style,
        key: container.getAttribute('data-tiles-key') || '',
        maxZoom: 19,
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors ' +
          '&copy; <a href="https://carto.com/attributions">CARTO</a>'
      });
    }
    return tileLayers[style];
  }

  var currentStyle = tileStyle();
  tileLayer(currentStyle).addTo(map);

  // Follow the site's light/dark toggle. Swap layers rather than calling redraw(),
  // which asks for fractional zoom tiles (light_all/1.5/...) when the map sits between
  // zoom levels. The new layer loads on top, then the old one goes, so nothing flashes.
  new MutationObserver(function() {
    var style = tileStyle();
    if (style === currentStyle) return;
    var previous = tileLayer(currentStyle);
    var next = tileLayer(style).addTo(map);
    currentStyle = style;
    next.once('load', function() {
      if (currentStyle === style) map.removeLayer(previous);
    });
  }).observe(root, { attributes: true, attributeFilter: ['data-theme'] });

  // Until the map is clicked, the scroll wheel scrolls the page; say so when someone tries.
  var hint = L.DomUtil.create('div', 'talks-map-hint', container);
  hint.setAttribute('aria-hidden', 'true');
  hint.textContent = 'Click the map to zoom with the scroll wheel';
  var hintTimer;

  function hideHint() {
    clearTimeout(hintTimer);
    hint.classList.remove('is-visible');
  }

  if (!touchFirst) {
    container.addEventListener('wheel', function() {
      if (map.scrollWheelZoom.enabled()) return;
      hint.classList.add('is-visible');
      clearTimeout(hintTimer);
      hintTimer = setTimeout(hideHint, 1600);
    }, { passive: true });
  }

  function setActive(active) {
    if (active) {
      hideHint();
      map.scrollWheelZoom.enable();
      map.dragging.enable();
    } else {
      map.scrollWheelZoom.disable();
      if (touchFirst) map.dragging.disable();
    }
  }
  container.addEventListener('focusin', function() { setActive(true); });
  container.addEventListener('focusout', function(event) {
    if (!container.contains(event.relatedTarget)) setActive(false);
  });

  // The opening view: every pin, kept on the tiles (they end at 85 degrees). This is
  // not a hard limit: a tall popup near the top may still pan the map past the edge.
  var homeBounds = places.map(function(place) { return [place.lat, place.lng]; });
  var TILE_BOUNDS = [[-85, -540], [85, 540]];

  function showAll(animate) {
    if (homeBounds.length) {
      map.fitBounds(homeBounds, { padding: [40, 40], maxZoom: 5, animate: animate });
    } else {
      map.setView([25, 0], 1, { animate: animate });
    }
    if (animate) {
      map.once('moveend', function() { map.panInsideBounds(TILE_BOUNDS, { animate: false }); });
    } else {
      map.panInsideBounds(TILE_BOUNDS, { animate: false });
    }
  }

  // Leaflet only draws markers once the map has a view, so set it first.
  showAll(false);

  // "Show all pins" button: under the zoom buttons on desktop, bottom left on phones.
  var FRAME_ICON = '<svg viewBox="0 0 16 16" width="14" height="14" aria-hidden="true" focusable="false">' +
    '<path d="M2.5 6V2.5H6M10 2.5h3.5V6M13.5 10v3.5H10M6 13.5H2.5V10" fill="none" stroke="currentColor" ' +
    'stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/></svg>';

  function addResetButton(parent) {
    var button = L.DomUtil.create('a', 'talks-map-reset', parent);
    button.href = '#';
    button.title = 'Show all pins';
    button.setAttribute('role', 'button');
    button.setAttribute('aria-label', 'Show all pins');
    button.innerHTML = FRAME_ICON;
    L.DomEvent.disableClickPropagation(button);
    L.DomEvent.on(button, 'click', function(event) {
      L.DomEvent.preventDefault(event);
      map.closePopup();
      showAll(!reduceMotion);
    });
  }

  if (map.zoomControl) {
    addResetButton(map.zoomControl.getContainer());
  } else {
    var ResetControl = L.Control.extend({
      options: { position: 'bottomleft' },
      onAdd: function() {
        var bar = L.DomUtil.create('div', 'leaflet-bar');
        addResetButton(bar);
        return bar;
      }
    });
    new ResetControl().addTo(map);
  }

  function cityName(city) {
    return String(city).split(',')[0];
  }

  // Nearby cities merge into one numbered pin until the map is zoomed in far enough to
  // tell them apart (Leaflet.markercluster). Without the plugin, every pin shows on its own.
  var clusters = typeof L.markerClusterGroup === 'function';
  var pinLayer = clusters ? L.markerClusterGroup({
    // Wider when zoomed out, so the US East Coast folds into one group, not two touching ones.
    maxClusterRadius: function(zoom) { return zoom <= 3 ? 56 : 36; },
    showCoverageOnHover: false,
    spiderfyOnMaxZoom: false,
    zoomToBoundsOnClick: false,
    animate: !reduceMotion,
    iconCreateFunction: function(cluster) {
      var total = 0;
      var hasTalk = false;
      var names = [];
      cluster.getAllChildMarkers().forEach(function(marker) {
        total += marker.options.eventCount;
        hasTalk = hasTalk || marker.options.hasTalk;
        names.push(cityName(marker.options.city));
      });
      var className = 'talks-map-pin talks-map-pin--count talks-map-pin--cluster' + (hasTalk ? '' : ' talks-map-pin--attended');
      return L.divIcon({
        className: 'talks-map-marker talks-map-marker--cluster',
        html: '<span class="' + className + '" aria-hidden="true">' + total + '</span>' +
          '<span class="visually-hidden">' + escapeHtml(total + ' events in ' + names.join(', ') + '. Select to zoom in.') + '</span>',
        iconSize: [36, 36],
        iconAnchor: [18, 18],
        tooltipAnchor: [0, -14]
      });
    }
  }) : L.layerGroup();
  pinLayer.addTo(map);

  if (clusters) {
    pinLayer.on('clusterclick', function(event) {
      event.layer.zoomToBounds({ padding: [60, 60] });
    });
    if (!touchFirst) {
      // Name the cities in a group on hover, like a single pin names its city.
      pinLayer.on('clustermouseover', function(event) {
        var cluster = event.layer;
        if (!cluster.getTooltip()) {
          var names = cluster.getAllChildMarkers().map(function(marker) { return cityName(marker.options.city); });
          cluster.bindTooltip(names.join(' · '), { direction: 'top', className: 'talks-map-tooltip', opacity: 1 });
        }
        cluster.openTooltip();
      });
    }
  }

  // Popup size follows the map's size, so it is recomputed when the map resizes
  // (a phone turned sideways, a window made narrower).
  function popupWidth() { return Math.max(200, Math.min(300, container.clientWidth - 64)); }
  function popupMaxHeight() { return Math.max(200, container.clientHeight - 48); }
  // Keep popups clear of the zoom buttons (top left) and the attribution (bottom right).
  var panTopLeft = map.zoomControl ? [60, 12] : [12, 12];
  var panBottomRight = [12, 28];

  var markers = places.map(function(place) {
    var count = place.events.length;
    var size = count > 1 ? 30 : 24; // hit area; the visible pin is smaller
    var hasTalk = place.events.some(function(event) { return !event.attended; });
    var label = place.city + ': ' + count + (count > 1 ? ' events' : ' event');
    var marker = L.marker([place.lat, place.lng], {
      riseOnHover: true,
      // Talk pins sit above attended-only pins where they overlap.
      zIndexOffset: hasTalk ? 1000 : 0,
      city: place.city,
      eventCount: count,
      hasTalk: hasTalk,
      icon: L.divIcon({
        className: 'talks-map-marker',
        html: pinHtml(place),
        iconSize: [size, size],
        iconAnchor: [size / 2, size / 2],
        popupAnchor: [0, -size / 2 + 4],
        tooltipAnchor: [0, -size / 2 + 4]
      })
    });
    // The pin's element exists only while it is on the map (not folded into a group).
    marker.on('add', function() {
      marker.getElement().setAttribute('aria-label', label);
    });
    if (!touchFirst) {
      marker.bindTooltip(place.city, { direction: 'top', className: 'talks-map-tooltip', opacity: 1 });
    }
    var popup = buildPopup(place, function() { return marker.getPopup(); });
    marker.on('popupopen', function() {
      marker.closeTooltip();
      var element = marker.getElement();
      if (element) element.classList.add('is-active');
    });
    marker.on('popupclose', function() {
      var element = marker.getElement();
      if (element) element.classList.remove('is-active');
      popup.reset(); // so the next open starts on the newest event
    });
    marker.bindPopup(popup.element, {
      minWidth: popupWidth(),
      maxWidth: popupWidth(),
      maxHeight: popupMaxHeight(),
      autoPanPaddingTopLeft: panTopLeft,
      autoPanPaddingBottomRight: panBottomRight
    });
    return marker;
  });

  map.on('resize', function() {
    var width = popupWidth();
    var maxHeight = popupMaxHeight();
    markers.forEach(function(marker) {
      var popup = marker.getPopup();
      popup.options.minWidth = popup.options.maxWidth = width;
      popup.options.maxHeight = maxHeight;
      if (popup.isOpen()) popup.update();
    });
  });

  if (clusters) {
    pinLayer.addLayers(markers);
  } else {
    markers.forEach(function(marker) { pinLayer.addLayer(marker); });
  }
})();
