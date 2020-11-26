var gisaid = gisaid || (function () {
	
	this.logging  = false;
	this.path 	  = "js/gisaid"
	this.ackCache = {};

	var log = function(txt) {
		if(this.logging) console.log("[gisaid] " + txt);
	}
	
	var appendCssLink = function() {
		var head  = document.getElementsByTagName('head')[0];
		var body  = document.getElementsByTagName('body')[0];
		var link  = document.createElement('link');

		link.rel  = 'stylesheet';
		link.type = 'text/css';
		link.href = this.path + '/gisaid.css';

		if(head) {
			head.appendChild(link);
		} else {
			body.appendChild(link);
		}
	}
	
	var createPopupObject = function() {
		var body = document.getElementsByTagName('body')[0];
		fetch(this.path + '/popup.html')
		  .then(response => response.text())
		  .then(popupInnerHTML => { 
				window.popupDivElement = document.createElement("div");
				window.popupDivElement.classList.add("gisaid-popup");
				window.popupDivElement.innerHTML = popupInnerHTML;
				window.popupDivElement.hide = function() {
					this.style.visibility = "hidden";
					this.style.pointerEvents = "none";
				}
				window.popupDivElement.show = function() {
					this.style.visibility = "visible";
					this.style.pointerEvents = "all";
				}
			});
	}
	
	var generateAckLinkFromAccessionId = function(accessionId) {
		let accessionIdNumber = accessionId.split("_")[2];
		
		return "https://www.epicov.org/acknowledgement/" + accessionIdNumber.charAt(2)
				+ accessionIdNumber.charAt(3) + "/" + accessionIdNumber.charAt(4)
				+ accessionIdNumber.charAt(5) + "/" + accessionId + ".json"
	}
	
	var hidePopup = function(wrappedElement) {
		let popupContainer = wrappedElement.parentNode.querySelector(".gisaid-popup");
		if(popupContainer) {
			console.log("[log] Hiding popup");
			popupContainer.style.pointerEvents = 'none';
			popupContainer.style.visibility = 'hidden';
		} else {
			setTimeout(function() {
				hidePopup(wrappedElement);				
			},120);
		}		
	}

	var renderPopup = function(wrappedElement) {

		let popupInstance = wrappedElement.querySelector(".gisaid-popup");
		let accessionId = wrappedElement.textContent.trim();
		if(wrappedElement.firstChild && wrappedElement.firstChild.getAttribute && wrappedElement.firstChild.getAttribute("epi_isl_id")) {
			accessionId = wrappedElement.firstChild.getAttribute("epi_isl_id");
		}
		let ackLink = generateAckLinkFromAccessionId(accessionId);	
	
		function positionPopup() {
			
			function getCurrentScrollPosition() {
				return {"y":document.documentElement.scrollTop,"x":document.documentElement.scrollLeft};	
			}

			let wrappedElementBoundingRectangle 	= wrappedElement.getBoundingClientRect();
			let popupBoundingRectangle 				= popupInstance.getBoundingClientRect();
			let windowWidth 			 			= window.innerWidth;
			let windowHeight 			 			= window.innerHeight;	
			let documentHorizontalCenter 			= (windowWidth/2);

			if(wrappedElementBoundingRectangle.x > (windowWidth/2)) {
				popupInstance.style.marginLeft = "-" + (popupBoundingRectangle.width - wrappedElementBoundingRectangle.width - 1) + "px";
			} else {
				popupInstance.style.marginLeft = "-1px";
			}
			if(wrappedElementBoundingRectangle.y > (windowHeight/2)) {
				popupInstance.style.marginTop = "-" + (popupBoundingRectangle.height) + "px";
			} else {
				popupInstance.style.marginTop = (wrappedElementBoundingRectangle.height) + "px";
			}		
			
			popupInstance.style.visibility = "visible";
			popupInstance.style.pointerEvents = "all";
			
		}
	
		if(!this.ackCache[ackLink]) {
			// fetch ack data and cache it		
			fetch(ackLink)
			  .then(response => response.json())
			  .then(ackData => {
				ackData["covv_accession_id"] = accessionId;   
				this.ackCache[ackLink] = ackData ;
				popupInstance = window.popupDivElement.cloneNode(true);
				popupInstance.style.fontSize = "inherit";
				popupInstance.accessionId = accessionId;		
				Object.keys(ackData).forEach(function(key,i) {
					popupInstance.innerHTML = popupInstance.innerHTML.replace("{" + key + "}",ackData[key]);
				});
				popupInstance.addEventListener("touchend",function(event) {
					hidePopup(event.target.parentNode.querySelector("div.gisaid-sequence"));		
				},true);

				if(!wrappedElement.parentNode.querySelector(".gisaid-popup")) {
					console.log("[log] Fetch: No popup there!")
					wrappedElement.parentNode.insertBefore(popupInstance, wrappedElement.parentNode.firstChild);
				} else {
					console.log("[log] Fetch: popup already there!");
				}
				
				positionPopup();
			  });
		} else {
			// get cached ack data
			var ackData = this.ackCache[ackLink];
			//(wrappedElement.parentNode.firstChild.classList + "").indexOf("gisaid-popup") == -1
			if(!wrappedElement.parentNode.querySelector(".gisaid-popup")) {
				console.log("[log] Cache popup not there!");
				popupInstance = window.popupDivElement.cloneNode(true);
				popupInstance.style.fontSize = "inherit";			
				Object.keys(ackData).forEach(function(key,i) {
					popupInstance.innerHTML = popupInstance.innerHTML.replace("{" + key + "}",ackData[key]);
				});
				popupInstance.addEventListener("touchend",function(event) {
					hidePopup(event.target.parentNode.querySelector("div.gisaid-sequence"));		
				},true);
				wrappedElement.parentNode.insertBefore(popupInstance, wrappedElement.parentNode.firstChild);
			} else {
				console.log("[log] Cache popup already there!");
				popupInstance = wrappedElement.parentNode.firstChild;
			}
			positionPopup();						  
		}			
	}
	
	var getAckData = function(accessionId,cb) {
		if(accessionId.indexOf("EPI_ISL_") === 0) {
			let ackLink = generateAckLinkFromAccessionId(accessionId);		
			if(!this.ackCache[ackLink]) {	
				fetch(ackLink)
				  .then(response => response.json())
				  .then(ackData => {
					  this.ackCache[ackLink] = ackData ;
					  cb(ackData);
				});
			} else {
				cb(this.ackCache[ackLink]);
			}
		} else {
			cb({});
		}	
	}
	 	
   	var addPopups = function() {
	   	
	   	[].forEach.call(document.querySelectorAll('body *[epi_isl_id]'),
					function(currentNode) {

	   		if (currentNode.getAttribute("epi_isl_id") !== "") {
	   			var gisaidContainer = document.createElement("div");
					gisaidContainer.classList.add("gisaid-container");

					var morphStyle = window.getComputedStyle(currentNode, null)
																 .getPropertyValue("display");

					if (morphStyle.indexOf("inline") > -1) {
						gisaidContainer.style.display = "inline-block";
					}
					else if(morphStyle.indexOf("block") === 0) {
						gisaidContainer.style.display = "block";
					}

					var gisaidSequence  = document.createElement("div");

					gisaidSequence.classList.add("gisaid-sequence");
					gisaidSequence.addEventListener("mouseover",function(e) {
						gisaid.renderPopup(this)
					} );
					gisaidSequence.addEventListener("mouseout",function(e)  {
						gisaid.hidePopup(this)
					} );

					var clonedNode = currentNode.cloneNode(true);

					gisaidSequence.appendChild(clonedNode);
					gisaidContainer.appendChild(gisaidSequence);

					currentNode.parentNode.replaceChild(gisaidContainer, currentNode);
				}
		}); 	   	
	}
	
	return {
		init:function() {
			appendCssLink();
			createPopupObject();
		},
		renderPopup:function(wrappedElement) {
			renderPopup(wrappedElement);
		},
		hidePopup:function(wrappedElement) {
			hidePopup(wrappedElement);		
		},
		addPopups:function() {
			addPopups();
		},
		getAcknowledgementData:function(epiIslId, cb) {
			getAckData(epiIslId,function(ackData) {
				cb(ackData);
			});
		}	
	}
})();

window.addEventListener('DOMContentLoaded', gisaid.init, false);
