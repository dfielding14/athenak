// Custom JavaScript for AthenaK Documentation

// Function to get the current theme
function getCurrentTheme() {
    const theme = document.body.dataset.theme;
    if (theme === 'dark' || 
        (theme === 'auto' && window.matchMedia('(prefers-color-scheme: dark)').matches)) {
        return 'dark';
    }
    return 'light';
}

// Function to get theme-appropriate colors
function getMermaidTheme() {
    const isDark = getCurrentTheme() === 'dark';
    
    if (isDark) {
        return {
            theme: 'dark',
            themeVariables: {
                primaryColor: '#667eea',
                primaryTextColor: '#e8e6e3',
                primaryBorderColor: '#8b7ff7',
                lineColor: '#a084ca',
                secondaryColor: '#2d2d2d',
                tertiaryColor: '#1a1a1a',
                background: '#1e1e1e',
                mainBkg: '#667eea',
                secondBkg: '#2d2d2d',
                tertiaryBkg: '#1a1a1a',
                textColor: '#e8e6e3',
                labelTextColor: '#e8e6e3',
                nodeTextColor: '#e8e6e3',
                edgeLabelBackground: '#2d2d2d',
                clusterBkg: '#2d2d2d',
                clusterBorder: '#4a4a4a',
                defaultLinkColor: '#a084ca',
                titleColor: '#e8e6e3',
                actorTextColor: '#1a1a1a',
                actorLineColor: '#8b7ff7'
            }
        };
    } else {
        return {
            theme: 'neutral',
            themeVariables: {
                primaryColor: '#667eea',
                primaryTextColor: '#333',
                primaryBorderColor: '#667eea',
                lineColor: '#764ba2',
                secondaryColor: '#f0f0f0',
                tertiaryColor: '#fff',
                background: '#ffffff',
                mainBkg: '#667eea',
                secondBkg: '#f0f0f0',
                tertiaryBkg: '#ffffff',
                textColor: '#333',
                labelTextColor: '#333',
                nodeTextColor: '#333',
                edgeLabelBackground: '#f0f0f0',
                clusterBkg: '#f9f9f9',
                clusterBorder: '#ddd',
                defaultLinkColor: '#764ba2',
                titleColor: '#333'
            }
        };
    }
}

// Function to reinitialize Mermaid with current theme
function reinitializeMermaid() {
    if (typeof mermaid !== 'undefined') {
        const themeConfig = getMermaidTheme();
        
        // Initialize mermaid with theme-appropriate config
        mermaid.initialize({
            startOnLoad: false,
            theme: themeConfig.theme,
            themeVariables: themeConfig.themeVariables,
            flowchart: {
                useMaxWidth: true,
                htmlLabels: true,
                curve: 'basis'
            }
        });
        
        // Re-render all mermaid blocks
        const mermaidBlocks = document.querySelectorAll('.mermaid');
        if (mermaidBlocks.length > 0) {
            // Clear existing SVGs
            mermaidBlocks.forEach(function(block) {
                const svg = block.querySelector('svg');
                if (svg) {
                    svg.remove();
                }
                // Reset the block
                block.removeAttribute('data-processed');
            });
            
            // Re-render with new theme
            mermaid.init(undefined, mermaidBlocks);
        }
    }
}

// Enhanced collapsible sidebar function
function makeCollapsible(caption, list, isMainContent = false) {
    console.log(`makeCollapsible called:`, {
        isMainContent,
        captionHTML: caption.outerHTML.substring(0, 100),
        alreadyProcessed: caption.hasAttribute('data-collapsible-processed')
    });
    
    // Check which type of processing
    const processedAttr = isMainContent ? 'data-main-collapsible-processed' : 'data-sidebar-collapsible-processed';
    
    // Prevent duplicate processing
    if (caption.hasAttribute(processedAttr)) {
        console.log('  - Already processed, skipping');
        return false;
    }
    caption.setAttribute(processedAttr, 'true');
    
    const captionText = caption.querySelector('.caption-text');
    if (!captionText) {
        console.log('  - No caption-text found');
        return false;
    }
    
    // Style the caption as clickable
    caption.style.cursor = 'pointer';
    caption.style.userSelect = 'none';
    caption.style.position = 'relative';
    
    // Create arrow element
    const arrow = document.createElement('span');
    arrow.className = isMainContent ? 'main-caption-arrow' : 'caption-arrow';
    arrow.style.cssText = `
        display: inline-block;
        margin-right: 0.5em;
        transition: transform 0.2s;
        font-size: 0.8em;
        width: 0.8em;
        font-family: monospace;
    `;
    arrow.textContent = '▶';
    
    // Insert arrow before caption text
    captionText.insertBefore(arrow, captionText.firstChild);
    
    // Determine initial state
    let isExpanded;
    if (isMainContent) {
        // Main content starts collapsed
        isExpanded = false;
    } else {
        // Sidebar: expand if contains current page
        isExpanded = list.querySelector('.current') !== null;
    }
    
    // Set initial display state
    list.style.display = isExpanded ? 'block' : 'none';
    arrow.style.transform = isExpanded ? 'rotate(90deg)' : 'rotate(0deg)';
    caption.setAttribute('data-collapsed', !isExpanded ? 'true' : 'false');
    
    // Add click handler
    const toggleHandler = function(e) {
        e.preventDefault();
        e.stopPropagation();
        
        const isCollapsed = caption.getAttribute('data-collapsed') === 'true';
        
        if (isCollapsed) {
            list.style.display = 'block';
            arrow.style.transform = 'rotate(90deg)';
            caption.setAttribute('data-collapsed', 'false');
        } else {
            list.style.display = 'none';
            arrow.style.transform = 'rotate(0deg)';
            caption.setAttribute('data-collapsed', 'true');
        }
    };
    
    caption.addEventListener('click', toggleHandler);
    
    // Add hover effect
    caption.addEventListener('mouseenter', function() {
        arrow.style.color = 'var(--color-brand-primary, #667eea)';
    });
    
    caption.addEventListener('mouseleave', function() {
        arrow.style.color = '';
    });
    
    return true;
}

// Initialize collapsible sidebar
function initCollapsibleSidebar() {
    console.log('Initializing collapsible sidebar...');
    
    // Find sidebar container
    const sidebarContainers = [
        '.sidebar-tree',
        '.sidebar-scroll',
        '.sidebar-container',
        '[role="navigation"]'
    ];
    
    let sidebar = null;
    for (const selector of sidebarContainers) {
        sidebar = document.querySelector(selector);
        if (sidebar) {
            console.log('Found sidebar with selector:', selector);
            break;
        }
    }
    
    if (!sidebar) {
        console.log('No sidebar found');
        return;
    }
    
    // Process all captions in the sidebar
    const captions = sidebar.querySelectorAll('p.caption[role="heading"]');
    console.log('Found', captions.length, 'caption groups');
    
    let successCount = 0;
    captions.forEach((caption) => {
        // Find the associated list
        let list = caption.nextElementSibling;
        while (list && list.tagName !== 'UL') {
            list = list.nextElementSibling;
        }
        
        if (list) {
            if (makeCollapsible(caption, list, false)) {
                successCount++;
                const text = caption.querySelector('.caption-text')?.textContent || 'unknown';
                console.log('✓ Made collapsible:', text);
            }
        }
    });
    
    console.log('Successfully processed', successCount, 'of', captions.length, 'captions');
    
    // Also handle second-level items
    initSecondLevelCollapsible(sidebar);
}

// Initialize second-level collapsible items
function initSecondLevelCollapsible(container) {
    const l1Items = container.querySelectorAll('li.toctree-l1');
    console.log('Processing', l1Items.length, 'L1 items for second-level collapsibility');
    
    l1Items.forEach((item) => {
        const nestedList = item.querySelector('ul');
        if (!nestedList || nestedList.children.length === 0) {
            return;
        }
        
        const link = item.querySelector(':scope > a');
        if (!link || link.querySelector('.l2-arrow')) {
            return;
        }
        
        // Add arrow
        const arrow = document.createElement('span');
        arrow.className = 'l2-arrow';
        arrow.style.cssText = `
            display: inline-block;
            margin-right: 0.5em;
            transition: transform 0.2s;
            font-size: 0.8em;
            cursor: pointer;
        `;
        arrow.textContent = '▶';
        
        link.insertBefore(arrow, link.firstChild);
        
        // Set initial state
        // In sidebar, expand if current; in main content, start collapsed
        const inSidebar = container.classList.contains('sidebar-tree') || container.closest('.sidebar-tree');
        const hasCurrentChild = nestedList.querySelector('.current');
        const isExpanded = inSidebar ? (hasCurrentChild || item.classList.contains('current')) : false;
        
        nestedList.style.display = isExpanded ? 'block' : 'none';
        arrow.style.transform = isExpanded ? 'rotate(90deg)' : 'rotate(0deg)';
        item.setAttribute('data-l2-collapsed', !isExpanded ? 'true' : 'false');
        
        // Click handler for arrow only
        arrow.addEventListener('click', function(e) {
            e.preventDefault();
            e.stopPropagation();
            
            const isCollapsed = item.getAttribute('data-l2-collapsed') === 'true';
            
            if (isCollapsed) {
                nestedList.style.display = 'block';
                arrow.style.transform = 'rotate(90deg)';
                item.setAttribute('data-l2-collapsed', 'false');
            } else {
                nestedList.style.display = 'none';
                arrow.style.transform = 'rotate(0deg)';
                item.setAttribute('data-l2-collapsed', 'true');
            }
        });
        
        // Hover effect
        arrow.addEventListener('mouseenter', function() {
            arrow.style.color = 'var(--color-brand-primary, #667eea)';
        });
        
        arrow.addEventListener('mouseleave', function() {
            arrow.style.color = '';
        });
    });
}

// Initialize main content collapsible sections
function initCollapsibleMainContent() {
    console.log('Initializing main content collapsible sections...');
    
    // Try multiple selectors for main content area
    const contentSelectors = ['.content', 'main', 'article', '[role="main"]', '.document'];
    let mainContent = null;
    
    for (const selector of contentSelectors) {
        mainContent = document.querySelector(selector);
        if (mainContent) {
            console.log('Found main content with selector:', selector);
            break;
        }
    }
    
    if (!mainContent) {
        console.log('ERROR: No main content area found!');
        return;
    }
    
    // Find all toctree wrappers
    const wrappers = mainContent.querySelectorAll('.toctree-wrapper.compound');
    console.log('Found', wrappers.length, 'toctree sections in main content');
    
    let processedCount = 0;
    wrappers.forEach((wrapper, index) => {
        console.log(`Processing wrapper ${index}:`);
        const caption = wrapper.querySelector('p.caption[role="heading"]');
        const list = wrapper.querySelector('ul');
        
        if (caption && list) {
            console.log('  - Has caption and list');
            if (makeCollapsible(caption, list, true)) {
                processedCount++;
                const text = caption.querySelector('.caption-text')?.textContent || 'unknown';
                console.log('  ✓ Made main content collapsible:', text);
            } else {
                console.log('  ✗ Failed to make collapsible');
            }
        } else {
            console.log('  - Missing caption or list:', {hasCaption: !!caption, hasList: !!list});
        }
        
        // Also handle second-level in main content
        if (list) {
            initSecondLevelCollapsible(wrapper);
        }
    });
    
    console.log(`Processed ${processedCount} of ${wrappers.length} main content sections`);
}

// Robust initialization with retry mechanism
function initializeCollapsibles() {
    console.log('Starting collapsible initialization...');
    
    // Try sidebar
    initCollapsibleSidebar();
    
    // Try main content
    initCollapsibleMainContent();
    
    // Schedule a retry if needed
    setTimeout(() => {
        console.log('Checking for unprocessed elements...');
        
        // Check sidebar
        const sidebarCaptions = document.querySelectorAll('.sidebar-tree p.caption[role="heading"]:not([data-collapsible-processed])');
        if (sidebarCaptions.length > 0) {
            console.log('Found', sidebarCaptions.length, 'unprocessed sidebar captions, retrying...');
            initCollapsibleSidebar();
        }
        
        // Check main content
        const mainCaptions = document.querySelectorAll('.content p.caption[role="heading"]:not([data-collapsible-processed])');
        if (mainCaptions.length > 0) {
            console.log('Found', mainCaptions.length, 'unprocessed main captions, retrying...');
            initCollapsibleMainContent();
        }
    }, 500);
}

// Main initialization
document.addEventListener('DOMContentLoaded', function() {
    console.log('DOM loaded - initializing AthenaK documentation features');
    
    // Initialize collapsibles
    initializeCollapsibles();
    
    // Initialize Mermaid
    reinitializeMermaid();
    
    // Watch for theme changes
    const observer = new MutationObserver(function(mutations) {
        mutations.forEach(function(mutation) {
            if (mutation.type === 'attributes' && mutation.attributeName === 'data-theme') {
                setTimeout(reinitializeMermaid, 100);
            }
        });
    });
    
    observer.observe(document.body, {
        attributes: true,
        attributeFilter: ['data-theme']
    });
    
    // Listen for system theme changes
    if (window.matchMedia) {
        window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', function() {
            if (document.body.dataset.theme === 'auto') {
                setTimeout(reinitializeMermaid, 100);
            }
        });
    }
    
    // Additional retry after longer delay for slow-loading elements
    setTimeout(() => {
        console.log('Final initialization check...');
        initializeCollapsibles();
    }, 2000);
});

// Also try initialization when window fully loads (includes all resources)
window.addEventListener('load', function() {
    console.log('Window fully loaded - final initialization attempt');
    setTimeout(initializeCollapsibles, 100);
});