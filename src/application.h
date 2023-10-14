#pragma once
#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>
#include <vulkan/vulkan.h>
#include <vector>
#include <iostream>
#include <optional>

struct QueueFamilyIndices {
    std::optional<uint32_t> graphicsFamily;

    bool isComplete(){
        return graphicsFamily.has_value();
    }
};

class Application {
public:
    void run(int width, int height);

private:
    void initWindow();

    void checkAvailableExtensions();
    bool checkValidationLayerSupport();
    std::vector<const char*> getRequiredExtensions();
    void populateDebugMessengerCreateInfo(VkDebugUtilsMessengerCreateInfoEXT& createInfo);
    void setupDebugMessenger();

    void checkAvailableDevices(const std::vector<VkPhysicalDevice> &devices);
    bool isDeviceSuitable(VkPhysicalDevice device);
    void pickPhysicalDevice();

    QueueFamilyIndices findQueueFamilies(VkPhysicalDevice device);

    void createLogicalDevice();

    void createInstance();
    void initVulkan();

    void mainLoop();

    void cleanup();

    static VKAPI_ATTR VkBool32 VKAPI_CALL debugCallback(
        VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity,
        VkDebugUtilsMessageTypeFlagsEXT messageType,
        const VkDebugUtilsMessengerCallbackDataEXT* pCallbackData,
        void* pUserData) {

        std::cerr << "Validation layer: " << pCallbackData->pMessage << std::endl;

        return VK_FALSE;
    }

    GLFWwindow* m_window;

    VkInstance m_instance;
    VkPhysicalDevice m_physicalDevice = VK_NULL_HANDLE;
    VkPhysicalDeviceFeatures m_deviceFeatures{};
    VkDevice m_device = VK_NULL_HANDLE;
    VkQueue m_graphicsQueue;

    VkDebugUtilsMessengerEXT m_debugMessenger;

    int m_width, m_height;

    const std::vector<const char*> validationLayers = {
        "VK_LAYER_KHRONOS_validation"
    };

#ifdef NDEBUG
    const bool enableValidationLayers = false;
#else
    const bool enableValidationLayers = true;
#endif

};
